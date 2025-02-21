"""TinkerHap - Read-Based Phasing Algorithm with Integrated Multi-Method Support for Enhanced Accuracy"""

# Requirements: Python >=3.6.0, pysam >=0.17.0

from typing import Optional
import re
import time
import sys
import os
import argparse
import bisect
import pysam

APP_NAME = "TinkerHap"
APP_DESC = "Read-Based Phasing Algorithm with Integrated Multi-Method Support for Enhanced Accuracy"
APP_VERSION = "1.0.0"
APP_DATE = "2025/02/13"
APP_URL = "https://github.com/DZeevi-Lab/TinkerHap"

GRAVITY_SNP = 2
GRAVITY_INDEL = 1
OVERWRITE_OUTPUTS = True


def invert_phasing(phase: int) -> int:
    """Swaps the phase of an allele (1 <-> 2)."""
    return 1 if phase == 2 else 2


def phasing_by_linkage(phase: int, linkage: int) -> int:
    """Return phase or invert it based on linkage score."""
    return phase if linkage > 0 else invert_phasing(phase)


def fill_del(allele: str, ref: str) -> str:
    """Adjusts allele to match reference length, filling dels with dots."""
    if (allele == ref) and (len(ref) > 1):
        return allele[0]

    return allele.ljust(1+len(ref)-len(allele), ".")


def fill_dels(alleles: list, ref: str) -> list:
    """Formats both alleles based with fill_del on reference length."""
    return [fill_del(alleles[0], ref), fill_del(alleles[1], ref)]


def idx_from_positions(read_positions: list, position: int, start_index: int) -> tuple:
    """Finds index of a position in a sorted list, handling None values."""
    low = start_index
    high = len(read_positions) - 1

    while low <= high:
        mid = (low + high) // 2
        value = read_positions[mid]

        # Skip None values by adjusting the mid point
        if value is None:
            # Move to the next non-None element (linear step)
            left, right = mid - 1, mid + 1
            while left >= low and read_positions[left] is None:
                left -= 1
            while right <= high and read_positions[right] is None:
                right += 1

            if left >= low and read_positions[left] is not None:
                mid = left
                value = read_positions[mid]
            elif right <= high and read_positions[right] is not None:
                mid = right
                value = read_positions[mid]
            else:
                return (None, mid)  # No valid positions found

        if value == position:
            return (mid, mid)  # Found the desired position

        if value > position:
            high = mid - 1  # Search the left half
        else:
            low = mid + 1  # Search the right half

    return (None, low)  # If not found, return None and the next index to search from


class Read():
    """Lightweight structure for storing read data."""
    __slots__ = 'read', 'name', 'start', 'end', 'p1', 'p2', 'hp', 'ht', 'snps', 'alleles'


class BAMPhaser:
    """Handles read-based phasing using VCF and BAM data."""
    pair_end: bool = True
    sample_id: str = ""
    vcf_in_path: str = ""
    vcf_out_path: str = ""
    bam_in_path: str = ""
    bam_out_path: str = ""
    bed_in_path: str = ""
    bed_out_path: str = ""
    log_file_path: str = ""
    vcf_scaffold_path: str = ""
    input_region: str = ""
    quiet_mode: bool = True
    filter_map_quality: int = 20
    max_depth: int = 100

    log_file = None
    max_read_size: int
    progress_time: float = 0

    phase_count: int
    vcf_htz_sites: dict
    vcf_htz_sites_positions: list
    vcf_variants: list
    prev_phasing: tuple
    haplotype_start: int
    haplotype_end: int
    haplotype_read_idx: int
    haplotype_no: int
    vcf_writer: pysam.VariantFile
    vcf_reader: pysam.VariantFile
    bam_in: pysam.AlignmentFile
    bam_out: pysam.AlignmentFile
    bam_outx: pysam.AlignmentFile
    bam_out1: pysam.AlignmentFile
    bam_out2: pysam.AlignmentFile
    haplotypes: dict
    time_started: float
    haplotypes_scaffolds: dict
    log_file = None

    def print_progress(self, position: int, steps: int, comment: str = "") -> None:
        """Displays periodic progress updates to stderr."""
        if self.quiet_mode:
            return

        if ((position % steps == 0) or (time.time() - self.progress_time > 15)):
            self.progress_time = time.time()
            if comment != "":
                comment = " ("+str(comment)+")"
            print(f"  {position:,}{comment}", end="\r", file=sys.stderr)

    def log(self, text: str, print_out: bool = True) -> None:
        """Logs a message to a file and optionally prints it."""
        if not hasattr(self, "time_started"):
            self.time_started = time.time()
            if self.log_file_path != "":
                self.log_file = open(self.log_file_path, "w", encoding="utf-8")

        text = str(round(time.time() - self.time_started))+"] " + text

        if self.log_file is not None:
            self.log_file.write(text + "\n")
        if (print_out) and (not self.quiet_mode):
            print(text, file=sys.stderr)

    def parse_bed_file(self, bed_file: str) -> list:
        """Reads BED file and extracts chromosome regions."""
        bed_data = []
        with open(bed_file, 'r', encoding="utf-8") as file:
            for line in file:
                # Skip comments
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.strip().split()
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3]
                bed_data.append((chrom, start, end, name))
        return bed_data

    def calculate_read_phasing(self, read: Read) -> int:
        """Determines the haplotype phase of a read."""
        phase1 = read.p1
        phase2 = read.p2

        if phase1 > phase2:
            return 1
        if phase2 > phase1:
            return 2

        return 0

    def increase_phase_score(self, read: Read, phasing: int, gravity: int) -> None:
        """Updates a read's phase score with a given weight."""
        self.phase_count = max(self.phase_count, phasing)
        if phasing == 1:
            read.p1 += gravity
        else:
            read.p2 += gravity

    def link_reads(self, read1: Read, read2: Read, linkage: int, read1idx: int) -> None:
        """Links two reads by updating their haplotype phase scores."""

        gravity = abs(linkage)

        phasing1 = self.calculate_read_phasing(read1)
        phasing2 = self.calculate_read_phasing(read2)

        if (phasing1 == 0) and (phasing2 == 0):
            htz_sites_positions = list(read1.alleles)
            if (self.haplotype_start > 0) and (htz_sites_positions[0] > self.haplotype_end):
                self.new_haplotype(read1idx)

            phasing1 = 1
            phasing2 = phasing_by_linkage(phasing1, linkage)
        elif (phasing1 == 0):
            phasing1 = phasing_by_linkage(phasing2, linkage)
            phasing2 = phasing_by_linkage(phasing1, linkage)
        elif (phasing2 == 0):
            phasing2 = phasing_by_linkage(phasing1, linkage)
            phasing1 = phasing_by_linkage(phasing2, linkage)
        else:
            if (linkage > 0):
                if (phasing1 != phasing2):
                    if (phasing1 == 1) and (read1.p1 < read2.p2):
                        phasing1 = 2
                    elif (phasing1 == 2) and (read1.p2 < read2.p1):
                        phasing1 = 1
                phasing2 = phasing1
            else:
                if phasing1 == phasing2:
                    if (phasing1 == 1) and (read1.p1 < read2.p1):
                        phasing1 = 2
                        phasing2 = 1
                    elif (phasing1 == 2) and (read1.p2 < read2.p2):
                        phasing1 = 2
                        phasing2 = 1
                    else:
                        phasing2 = invert_phasing(phasing1)

                p2 = invert_phasing(phasing1)
                phasing1 = invert_phasing(phasing2)
                phasing2 = p2

        self.increase_phase_score(read1, phasing1, gravity)
        self.increase_phase_score(read2, phasing2, gravity)

    def get_allele(self, sequence: str, positions: list, idx: int) -> str:
        """Extracts the allele at a specific position in a read."""
        allele = sequence[idx]
        while idx+1 < len(positions):
            if positions[idx+1] is None:  # Insert
                allele += sequence[idx+1]
                idx += 1
            else:
                # When next position is more than +1 - it means we have a deletion
                if positions[idx] is not None:
                    allele += "."*(positions[idx+1]-positions[idx]-1)
                break

        return allele

    def htz_allele_exists(self, htz_site: tuple, sequence: str, positions: list, idx: int) -> str:
        """Checks if a heterozygous allele exists at a position."""
        allele = self.get_allele(sequence, positions, idx)
        if allele in (htz_site[0], htz_site[1]):
            return allele
        return "-"

    def read_name_reversed(self, read: Read) -> str:
        """Reverses the strand suffix in a read's name."""
        strand = read.name[-1]
        if strand == '+':
            strand = '-'
        else:
            strand = '+'

        return read.name[:-1]+strand

    def variant_is_snp(self, read: Read, position: int) -> bool:
        """Determines if a variant is an SNP or an indel."""
        return read.snps[position]

    def allele_exists(self, read: Read, site: int):
        """Checks if a specific allele exists in a read."""
        if site in read.alleles:
            return read.alleles[site]

        return False

    def test_htz_reads(self, read1: Read, read2: Read) -> int:
        """Computes linkage score between two reads."""
        if read1.start > read2.end:
            return 0

        if read1.end < read2.start:
            return 0

        reads_linked = 0

        for site in read1.alleles:
            allele1 = self.allele_exists(read1, site)
            if not allele1:
                continue

            allele2 = self.allele_exists(read2, site)
            if not allele2:
                continue

            gravity = GRAVITY_SNP if self.variant_is_snp(read1, site) else GRAVITY_INDEL

            reads_linked += gravity if allele1 == allele2 else -gravity

        return reads_linked

    def vcf_get_ref_end(self, contig: str) -> int:
        """Gets the end position of a reference contig from the VCF."""
        return self.vcf_reader.header.contigs[contig].length

    def vcf_read_region(self, sample_id: str, contig: str,
                        start: Optional[int], stop: Optional[int]) -> int:
        """Loads heterozygous sites from the VCF for a given region."""
        self.vcf_htz_sites = {}
        self.vcf_variants = []
        records_count = 0
        for record in self.vcf_reader.fetch(contig, start, stop):
            records_count += 1
            if records_count == 1:
                self.log(f"Reading variants for {contig}", True)

            self.vcf_variants.append(record)
            if record.qual == 0:
                continue

            if sample_id == "":
                data = record.samples[0]
            else:
                data = record.samples.get(sample_id)
            if data is None:
                continue

            # Check if site is HTZ
            if (data.alleles[0] != None) and (data.alleles[1] != None) and (data.alleles[0] != data.alleles[1]):
                alleles = data.alleles
                l1 = len(alleles[0])
                if l1 > 1:
                    alleles = (alleles[0][0]+("."*(l1-1)), alleles[1])  # DEL - fill with dots

                self.vcf_htz_sites[record.pos] = alleles

        self.vcf_htz_sites_positions = sorted(self.vcf_htz_sites.keys())
        return records_count

    def vcf_sample_rephase(self, sample: pysam.VariantRecordSample, phase: tuple) -> None:
        """Adjusts sample phasing in the VCF"""
        sample.alleles = phase
        sample.phased = True

    def vcf_count_alleles(self, sample: pysam.VariantRecordSample, reads: list,
                          variant: pysam.VariantRecord, reads_start: int,
                          allele_size: int = 1) -> tuple:
        """Counts allele occurrences across haplotypes for a given variant."""
        bam_pos = variant.pos - 1
        alleles_hp_count = [[0, 0], [0, 0]]  # Count of alleles per haplotype: alleles_hp_count[allele 0/1][haplotype]
        first_read = True
        haplotype = 0
        ref = variant.ref
        alleles = fill_dels(sample.alleles, ref)
        for idx in range(reads_start, len(reads)):
            read = reads[idx]
            if read.start - 1 > bam_pos:
                break

            if read.end+self.max_read_size < bam_pos:
                continue

            if first_read:
                reads_start = idx
                first_read = False

            if read.hp == 0:
                continue

            if not variant.pos in read.alleles:
                continue
            allele = read.alleles[variant.pos]

            haplotype = read.ht

            allele_no = -1
            if allele in alleles:
                allele_no = alleles.index(allele)
            elif allele_size == 1:
                if alleles[0][0] == ref[0]:
                    allele_no = 1
                else:
                    allele_no = 0

            if allele_no > -1:
                alleles_hp_count[allele_no][read.hp-1] += 1

        return (alleles_hp_count, haplotype, reads_start)

    def vcf_phase_htz_sample(self, sample: pysam.VariantRecordSample, reads: list,
                             variant: pysam.VariantRecord, reads_start: int) -> int:
        """Determines and records phasing information for a heterozygous site."""
        (alleles_hp_count, haplotype, reads_start) = self.vcf_count_alleles(sample, reads,
                                                                            variant, reads_start)

        if alleles_hp_count[0][1] == 0:
            ratio0 = alleles_hp_count[0][0]+1
        else:
            ratio0 = alleles_hp_count[0][0] / alleles_hp_count[0][1]

        if alleles_hp_count[1][1] == 0:
            ratio1 = alleles_hp_count[1][0]+1
        else:
            ratio1 = alleles_hp_count[1][0] / alleles_hp_count[1][1]

        allele0 = 0 if ratio0 > 1 else 1
        allele1 = 0 if ratio1 > 1 else 1

        if allele0 == allele1:
            if ratio0 > ratio1:
                allele0 = 0
                allele1 = 1
            elif (ratio0 == 0) and (ratio1 > 0):
                allele0 = 1
                allele1 = 0
            elif (((alleles_hp_count[0][0] == 0) and (alleles_hp_count[0][1] == 0)) or
                    ((alleles_hp_count[1][0] == 0) and (alleles_hp_count[1][1] == 0))):
                # When we cannot find an allele - assume it's the same phasing the previous one
                (allele0, allele1) = self.prev_phasing
            elif alleles_hp_count[0][0] > alleles_hp_count[1][0]:
                allele0 = 0
                allele1 = 1
            elif alleles_hp_count[0][0] < alleles_hp_count[1][0]:
                allele0 = 1
                allele1 = 0
            else:
                (allele0, allele1) = self.prev_phasing

        self.prev_phasing = (allele0, allele1)
        phase = (sample.alleles[allele0], sample.alleles[allele1])

        if allele0 != allele1:
            self.vcf_sample_rephase(sample, phase)
        else:
            sample.alleles = phase

        sample.phased = haplotype != 0
        sample['PS'] = haplotype if haplotype > 0 else None

        if haplotype in self.haplotypes_scaffolds:
            sample['PSS'] = self.haplotypes_scaffolds[haplotype]

        return reads_start

    def variant_get_sample(self, variant: pysam.VariantRecord, sample_id) -> pysam.VariantRecordSample:
        """Retrieves a sample from a VCF variant, using a sample ID if provided."""
        if sample_id == "":
            return variant.samples[0]
        else:
            return variant.samples.get(sample_id)

    def vcf_variant_revert(self, pending_variant: pysam.VariantRecord, pending_sample_data: list, sample_id: str) -> None:
        """Restores the original sample phasing in a VCF variant from stored data."""
        pending_sample = self.variant_get_sample(pending_variant, sample_id)
        pending_sample.alleles = pending_sample_data[0]
        pending_sample.phased = pending_sample_data[1]
        ps = pending_sample_data[2]
        if (ps is None):
            pending_sample['PS'] = None
        else:
            pending_sample['PS'] = ps

    def vcf_save_pending(self, pending_variant: pysam.VariantRecord, pending_undecided: list) -> None:
        """Writes the pending variant and undecided variants to the output VCF."""
        self.vcf_writer.write(pending_variant)
        for variant in pending_undecided:
            self.vcf_writer.write(variant)
        pending_undecided.clear()

    def vcf_save_reads(self, reads: list, sample_id: str) -> None:
        """Saves phased reads and variants back to the output VCF."""
        if not self.vcf_writer:
            return

        self.prev_phasing = (0, 1)
        self.log('Saving to VCF', True)

        current_haplotype = None
        pending_variant = None
        pending_undecided = []
        reads_start = 0
        last_pss = 0
        last_pss_haplotype = 0
        for variant in self.vcf_variants:
            sample = self.variant_get_sample(variant, sample_id)

            if (sample is None) or (sample.alleles[0] == None) or (sample.alleles[1] == None):
                if (pending_variant is None):
                    self.vcf_writer.write(variant)
                else:
                    pending_undecided.append(variant)
                continue

            org_sample_data = [sample.alleles, sample.phased, sample['PS'] if 'PS' in sample else None]
            if sample.alleles[0] == sample.alleles[1]:
                sample.phased = True
                sample['PS'] = current_haplotype
            else:
                reads_start = self.vcf_phase_htz_sample(sample, reads, variant, reads_start)

            if ('PS' in sample) and (sample['PS'] is not None):
                if ('PSS' in sample):
                    if last_pss != sample['PSS']:
                        last_pss = sample['PSS']
                        last_pss_haplotype = sample['PS']

                    org_sample_data[2] = last_pss_haplotype

                if (current_haplotype == sample['PS']):
                    if (pending_variant is not None):
                        self.vcf_save_pending(pending_variant, pending_undecided)
                        pending_undecided = []
                        pending_variant = None
                    self.vcf_writer.write(variant)
                else:
                    current_haplotype = sample['PS']
                    if (pending_variant is not None):
                        self.vcf_variant_revert(pending_variant, pending_sample_data, sample_id)
                        self.vcf_save_pending(pending_variant, pending_undecided)
                        pending_undecided = []
                    pending_sample_data = org_sample_data
                    pending_variant = variant
            else:
                if (pending_variant is None):
                    self.vcf_writer.write(variant)
                else:
                    pending_undecided.append(variant)

        if (pending_variant is not None):
            self.vcf_variant_revert(pending_variant, pending_sample_data, sample_id)
            self.vcf_save_pending(pending_variant, pending_undecided)

    def save_haplotypes(self) -> None:
        """Writes identified haplotype information to a BED file."""
        if self.bed_out_path == "":
            return

        with open(self.bed_out_path, "w", encoding="utf-8") as file:
            bedline = "\t".join(["#chrom", "start", "end", "haplotype"])
            file.write("#chrom\tstart\tend\thaplotype\n")
            for idx, haplotype in self.haplotypes.items():
                name = 'h-'+str(idx)
                bedline = "\t".join(["chr6", str(haplotype[0]), str(haplotype[1]), name])
                file.write(bedline+"\n")

    def htz_in_read(self, read_sequence: str, read_positions: list,
                    ref_start: int, ref_end: int) -> tuple:
        """Identifies heterozygous sites within a read's reference range."""
        alleles = {}
        snps = {}
        last_idx = 0

        htz_index = bisect.bisect_left(self.vcf_htz_sites_positions, ref_start)
        while htz_index < len(self.vcf_htz_sites_positions):
            position = self.vcf_htz_sites_positions[htz_index]
            if position > ref_end:
                break
            htz_index += 1

            (idx, last_idx) = idx_from_positions(read_positions, position-1, last_idx)
            if idx is None:
                continue

            htz_site = self.vcf_htz_sites[position]
            allele = self.htz_allele_exists(htz_site, read_sequence, read_positions, idx)
            if allele != "-":
                alleles[position] = allele
                snps[position] = (len(htz_site[0]) == 1) and (len(htz_site[1]) == 1)
        return (alleles, snps)

    def extend_haplotype(self, htz_sites: dict) -> None:
        """Expands the haplotype region based on heterozygous site bounds."""
        htz_sites_positions = list(htz_sites)
        if (self.haplotype_start == 0) or (self.haplotype_start > htz_sites_positions[0]):
            self.haplotype_start = htz_sites_positions[0]

        self.haplotype_end = max(self.haplotype_end, htz_sites_positions[-1])

    def phase_read_by_htz(self, reads: list, read1idx: int) -> None:
        """Assigns a phase to a read based on heterozygous site comparisons."""
        read = reads[read1idx]
        self.print_progress(read1idx, 1000, read.start)

        if len(read.alleles) > 0:
            reads_compared = 0
            for read2idx in range(read1idx+1, len(reads)):
                read2 = reads[read2idx]

                if read2.start > read.end:
                    break

                reads_compared += 1
                if reads_compared > self.max_depth:
                    break

                linkage = self.test_htz_reads(read, read2)

                if linkage != 0:
                    self.link_reads(read, read2, linkage, read1idx)

        phasing = self.calculate_read_phasing(read)
        read.hp = phasing
        read.ht = self.haplotype_no

        if phasing in (1, 2):
            self.extend_haplotype(read.alleles)
            if self.haplotype_read_idx == 0:
                self.haplotype_read_idx = read1idx
        elif (self.phase_count >= 2) and (len(read.alleles) > 0):
            self.new_haplotype(read1idx)

    def new_haplotype(self, readidx) -> None:
        """Starts a new haplotype and updates relevant metadata."""
        if self.haplotype_start > 0:
            self.haplotypes[self.haplotype_no] = [self.haplotype_start,
                                                  self.haplotype_end,
                                                  self.haplotype_read_idx]
            self.haplotype_no += 1
        self.phase_count = 0
        self.haplotype_start = 0
        self.haplotype_end = 0
        self.haplotype_read_idx = readidx

    def phase_by_htz(self, reads: list) -> None:
        """Performs initial phasing pass based on heterozygous site linkage."""
        self.log("Phasing by htz")

        self.haplotype_start = 0
        self.haplotype_end = 0
        self.haplotype_read_idx = 0
        for read1idx in range(len(reads)):
            self.phase_read_by_htz(reads, read1idx)

        if self.haplotype_start > 0:
            self.haplotypes[self.haplotype_no] = [self.haplotype_start,
                                                  self.haplotype_end,
                                                  self.haplotype_read_idx]
        self.log(f"HTZ: {len(self.haplotypes)} haplotypes")

    def haplotype_flip_direction(self, reads: list, ht: int) -> None:
        """Flips the phasing direction of a specified haplotype."""
        ht_found = False
        if ht not in self.haplotypes:
            return

        for idx in range(self.haplotypes[ht][2], len(reads)):
            if reads[idx].ht == ht:
                ht_found = True
                hp = reads[idx].hp
                if hp == 1:
                    hp = 2
                elif hp == 2:
                    hp = 1
            elif ht_found and (self.haplotypes[ht][1] < reads[idx].start+50):
                break

    def merge_haplotypes_direction(self, reads: list, ht1: int, ht2: int, direction: int) -> None:
        """Merges two haplotypes, flipping one if required by direction."""
        if (not ht1 in self.haplotypes) or (not ht2 in self.haplotypes):
            return

        ht_found = False
        for idx in range(self.haplotypes[ht2][2], len(reads)):
            if reads[idx].ht == ht2:
                ht_found = True
                reads[idx].ht = ht1
                if direction < 1:
                    hp = reads[idx].hp
                    if hp == 1:
                        hp = 2
                    elif hp == 2:
                        hp = 1
                    reads[idx].hp = hp
            elif ht_found and (self.haplotypes[ht2][1] < reads[idx].start+50):
                break

        if self.haplotypes[ht2][1] >= self.haplotypes[ht1][0]:
            self.haplotypes[ht1] = [self.haplotypes[ht1][0],
                                    self.haplotypes[ht2][1],
                                    self.haplotypes[ht1][2]]
            del self.haplotypes[ht2]
        else:
            raise ValueError(f"Unclear how to merge {ht1} and {ht2} at {self.haplotypes}")

    def merge_haplotypes(self, haplotypes_merging: dict, reads: list) -> None:
        """Merges haplotypes based on a merge map and adjusts read phasing."""
        self.log(f"Merging {len(haplotypes_merging)} / {len(self.haplotypes)}", True)

        for ht1 in haplotypes_merging:
            ht2 = ht1+1
            while ht2 in haplotypes_merging[ht1]:
                direction = haplotypes_merging[ht1][ht2]
                if direction == 0:
                    break
                direction_sign = 1 if direction > 0 else -1

                self.merge_haplotypes_direction(reads, ht1, ht2, direction)
                if ht2 in haplotypes_merging:
                    for key in haplotypes_merging[ht2]:
                        value = haplotypes_merging[ht2][key] * direction_sign
                        if key in haplotypes_merging[ht1]:
                            haplotypes_merging[ht1][key] += value
                        else:
                            haplotypes_merging[ht1][key] = value

                    haplotypes_merging[ht2] = {}
                ht2 += 1

    def phase_reads_pair(self, read1: Read, read2: Read, haplotypes_merging: dict) -> int:
        """Validates and updates haplotypes for paired-end reads."""
        pair_errors = 0

        read1_ht = read1.ht
        read1_hp = read1.hp
        read2_ht = read2.ht
        read2_hp = read2.hp

        if (read1_hp != read2_hp) and (read1_ht == read2_ht):
            if read1_hp == 0:
                read1.hp = read2_hp
            elif read2_hp == 0:
                read2.hp = read1_hp
            else:
                self.log(f"Pair error at {read1.name} ({read1.start})", False)
                pair_errors += 1

        if (read1_ht != read2_ht) and (read1_hp != 0) and (read2_hp != 0):
            if read1_ht > read2_ht:
                read1_ht = read2.ht
                read2_ht = read1.ht
            direction = 1 if read1_hp == read2_hp else -1

            if read1_ht in haplotypes_merging:
                if read2_ht in haplotypes_merging[read1_ht]:
                    haplotypes_merging[read1_ht][read2_ht] += direction
                else:
                    haplotypes_merging[read1_ht][read2_ht] = direction
            else:
                haplotypes_merging[read1_ht] = {read2_ht: direction}

        return pair_errors

    def phase_by_pair_ends(self, reads: list, read_names: dict) -> None:
        """Uses pair-end read information to refine haplotype phasing."""
        if not self.pair_end:
            return

        self.log('Phasing by pair-end', True)

        haplotypes_merging = {}
        pair_errors = 0

        for read1idx, read in enumerate(reads):
            self.print_progress(read1idx, 100)

            read1_name_r = self.read_name_reversed(read)
            if not read1_name_r in read_names:
                continue

            for read2idx in read_names[read1_name_r]:
                if read1_name_r == reads[read2idx].name:
                    pair_errors += self.phase_reads_pair(read, reads[read2idx], haplotypes_merging)

        self.merge_haplotypes(haplotypes_merging, reads)
        self.log(f"Total haplotypes: {len(self.haplotypes)}, Errors: {pair_errors}", True)

    def compare_alleles(self, allele1: list, allele2: list) -> int:
        """Compares two alleles and determines if they are identical or inverted."""
        if (allele1[0] == allele2[0]) and (allele1[1] == allele2[1]):
            return 1
        if (allele1[0] == allele2[1]) and (allele1[1] == allele2[0]):
            return -1
        return 0

    def scaffold_read_region(self, sample_id: str, contig: str, start: int, stop: int) -> dict:
        """Reads phased variants from a scaffold VCF for a specified region."""
        scaffold_sites = {}

        if not os.path.exists(self.vcf_scaffold_path):
            self.log("VCF file does not exist ("+self.vcf_scaffold_path+")")
            sys.exit(1)

        vcf_scaffold = pysam.VariantFile(self.vcf_scaffold_path)
        for record in vcf_scaffold.fetch(contig, start, stop):
            if sample_id == "":
                data = record.samples[0]
            else:
                data = record.samples.get(sample_id)

            if (data is None) or (not data.phased) or (data.alleles[0] == data.alleles[1]):
                continue

            alleles = data.alleles

            l1 = len(alleles[0])
            if l1 > 1:
                alleles = (alleles[0][0]+("."*(l1-1)), alleles[1])  # DEL - fill with dots

            scaffold_id = data['PS'] if 'PS' in data else 1
            scaffold_sites[record.pos] = (alleles, scaffold_id)

        vcf_scaffold.close()
        return scaffold_sites

    def phase_by_scaffold(self, reads: list, sample_id: str,
                          contig: str, start: int, stop: int) -> None:
        """Links haplotypes using phasing data from a scaffold VCF."""
        if self.vcf_scaffold_path == "":
            return

        if len(reads) == 0:
            return

        self.log('Phasing by scaffold', True)

        # HT might have changes for reads that has the same starting point as other - resort it
        reads.sort(key=lambda x: (x.ht, x.start))

        scaffold_sites = self.scaffold_read_region(sample_id, contig, start, stop)

        current_ht = reads[0].ht
        prev_ht = 0
        prev_ht_score = 0
        ht_score = 0
        scaffold_id = 0
        for read in reads:
            if current_ht != read.ht:
                if (ht_score != 0) and (scaffold_id != 0):
                    self.haplotypes_scaffolds[current_ht] = scaffold_id

                if ht_score < 0:
                    self.haplotype_flip_direction(reads, current_ht)

                if (prev_ht_score != 0) and (ht_score != 0):
                    self.haplotypes_scaffolds[current_ht] = scaffold_id
                    self.merge_haplotypes_direction(reads, prev_ht, current_ht, 1)
                else:
                    prev_ht = current_ht

                prev_ht_score = ht_score
                ht_score = 0
                scaffold_id = 0
                current_ht = read.ht

            for site in read.alleles:
                if site in scaffold_sites:
                    scaffold_id = max(scaffold_id, scaffold_sites[site][1])
                    scaffold_alleles = scaffold_sites[site][0]
                    if read.alleles[site] == scaffold_alleles[read.hp-1]:
                        score = 1
                    elif (read.alleles[site] == scaffold_alleles[0 if read.hp == 2 else 1]):
                        score = -1
                    else:
                        score = 0
                    ht_score += score

        self.log(f"Total haplotypes: {len(self.haplotypes)}", True)

    def load_reads(self, ref_contig: str, ref_start: int, ref_end: int) -> tuple:
        """Loads reads from the BAM file and processes their heterozygous sites."""
        self.max_read_size = 0
        reads = []
        read_names = {}
        current_pos = 0
        prev_loci = -1
        same_loci = 0
        for read in self.bam_in.fetch(ref_contig, ref_start, ref_end):
            self.print_progress(current_pos, 50000)

            if read.mapping_quality < self.filter_map_quality:
                continue

            positions = read.get_reference_positions(True)
            if len(positions) == 0:
                continue
            sequence = read.query_sequence

            if prev_loci == read.reference_start:
                same_loci += 1
                if same_loci > self.max_depth:
                    continue
            else:
                prev_loci = read.reference_start
                same_loci = 0

            htz_alleles, htz_snps = self.htz_in_read(sequence, positions,
                                                     read.reference_start, read.reference_end)
            if len(htz_alleles) == 0:
                continue

            allele_positions = list(htz_alleles)

            read_name = read.qname+('-' if read.is_reverse else '+')
            read_row = Read()
            read_row.name = read_name
            read_row.start = allele_positions[0]
            read_row.end = allele_positions[-1]
            read_row.p1 = 0
            read_row.p2 = 0
            read_row.hp = 0
            read_row.ht = 0

            read_row.alleles = htz_alleles
            read_row.snps = htz_snps

            if self.bam_out is not None:
                read_row.read = read

            self.max_read_size = max(self.max_read_size, read.reference_end - read.reference_start)
            reads.append(read_row)
            current_pos += 1

        reads.sort(key=lambda x: x.start)

        if self.pair_end:
            for current_pos, read in enumerate(reads):
                read_name = read.name
                if read_name in read_names:
                    read_names[read_name].append(current_pos)
                else:
                    read_names[read_name] = [current_pos]

        return (reads, read_names)

    def bam_save_reads(self, reads: list) -> None:
        """Saves phased reads to separate BAM files based on their haplotype."""
        if self.bam_out is None:
            return

        self.log('Saving to BAM', True)
        # Write reads to file
        for idx, read in enumerate(reads):
            self.print_progress(idx, 10000)

            read.read.set_tag('HP', read.hp)
            read.read.set_tag('HT', read.ht)

            self.bam_out.write(read.read)
            if read.hp == 1:
                self.bam_out1.write(read.read)
            elif read.hp == 2:
                self.bam_out2.write(read.read)
            else:
                self.bam_outx.write(read.read)

    def phase_region(self, sample_id: str, ref_contig: str,
                     ref_start: Optional[int], ref_end: Optional[int], region_name: str) -> None:
        """Phases reads in a specified genomic region."""
        self.haplotypes_scaffolds = {}
        self.phase_count = 0

        sample_hint = "\""+sample_id+"\"" if sample_id != "" else "<first>"

        if ref_start is None:
            ref_start = 0
        if ref_end is None:
            ref_end = self.vcf_get_ref_end(ref_contig)

        if self.vcf_read_region(sample_id, ref_contig, ref_start, ref_end + 189) == 0:
            return

        self.log(f"Loading reads of sample {sample_hint} for {region_name}")

        (reads, read_names) = self.load_reads(ref_contig, ref_start, ref_end)

        if len(self.vcf_htz_sites) > 0:
            reads_htz = f"{len(reads)/len(self.vcf_htz_sites):.2f}"
            bases_htz = f"{(ref_end-ref_start+1)/len(self.vcf_htz_sites):.2f}"
        else:
            reads_htz = "-"
            bases_htz = "-"

        self.log(f"{region_name} ({ref_contig}:{ref_start}-{ref_end}): {len(reads)} reads, " +
                 f"{len(self.vcf_htz_sites)} htz sites({reads_htz} reads/htz, " +
                 f"{bases_htz} bases/htz), largest read size: {self.max_read_size}")

        self.vcf_htz_sites = {}
        self.vcf_htz_sites_positions = []

        self.phase_by_htz(reads)
        self.phase_by_pair_ends(reads, read_names)
        self.bam_save_reads(reads)

        self.phase_by_scaffold(reads, sample_id, ref_contig, ref_start,
                               ref_end + self.max_read_size + 1)
        self.vcf_save_reads(reads, sample_id)

        self.log('Region done', True)

    def prepare_vcf_files(self):
        """Sets up VCF files for reading and writing phases."""
        if not os.path.exists(self.vcf_in_path):
            self.log("VCF file does not exist ("+self.vcf_in_path+")")
            sys.exit(1)

        self.vcf_reader = pysam.VariantFile(self.vcf_in_path)

        vcf_out_header = self.vcf_reader.header
        vcf_formats = vcf_out_header.formats.keys()
        if "PS" not in vcf_formats:
            vcf_out_header.formats.add("PS", 1, "Integer", "Phase set")
        if (self.vcf_scaffold_path != "") and ("PSS" not in vcf_formats):
            vcf_out_header.formats.add("PSS", 1, "Integer", "Phase set by scaffold")

        if self.vcf_out_path is not None:
            self.vcf_writer = pysam.VariantFile(self.vcf_out_path, mode="w", header=vcf_out_header)
        else:
            self.vcf_writer = None

    def prepare_bam_file(self, suffix: str) -> pysam.AlignmentFile:
        """Creates a new BAM file with the given suffix."""
        bam_dot = self.bam_out_path.rfind(".")
        return pysam.AlignmentFile(self.bam_out_path[:bam_dot] + suffix +
                                   self.bam_out_path[bam_dot:], "wb", template=self.bam_in)

    def prepare_bam_files(self):
        """Prepares BAM files for input and output operations."""

        if not os.path.exists(self.bam_in_path):
            self.log("BAM file does not exist ("+self.bam_in_path+")")
            sys.exit(1)

        self.bam_in = pysam.AlignmentFile(self.bam_in_path, "rb", require_index=True)
        if self.bam_out_path != "":
            self.bam_out = pysam.AlignmentFile(self.bam_out_path, "wb", template=self.bam_in)

            self.bam_outx = self.prepare_bam_file("X")
            self.bam_out1 = self.prepare_bam_file("1")
            self.bam_out2 = self.prepare_bam_file("2")
        else:
            self.bam_out = None

        self.bed_out_path = self.bed_out_path

    def phase_by_bed(self) -> None:
        """Phases reads using regions defined in an input BED file."""
        for bed_line in self.parse_bed_file(self.bed_in_path):
            self.phase_region(self.sample_id, bed_line[0], bed_line[1], bed_line[2], bed_line[3])

    def phase_all(self) -> None:
        """Phases all contigs in the VCF without region-specific constraints."""
        contigs = self.vcf_reader.header.contigs
        for contig in contigs:
            self.phase_region(self.sample_id, contig, None, None, "<all>")

    def phase_by_region(self) -> None:
        """Phases a user-defined region specified in the input arguments."""
        pattern = r"(?P<rname>[^:]+)(:(?P<startpos>\d+))?(?:-(?P<endpos>\d+))?"

        match = re.match(pattern, self.input_region)

        if match:
            rname = match.group('rname')
            startpos = int(match.group('startpos')) if match.group('startpos') else None
            endpos = int(match.group('endpos')) if match.group('endpos') else None

            self.phase_region(self.sample_id, rname, startpos, endpos, self.input_region)
        else:
            self.log("Invalid region format ("+self.input_region+")")
            sys.exit(1)

    def phase(self):
        """Phases data based on input files, regions, and configuration."""
        self.phase_count = 0
        self.haplotypes = {}
        self.haplotype_no = 1

        self.prepare_vcf_files()
        self.prepare_bam_files()

        if self.bed_in_path != "":
            self.phase_by_bed()
        elif self.input_region != "":
            self.phase_by_region()
        else:
            self.phase_all()

        self.save_haplotypes()

        if self.log_file is not None:
            self.log_file.close()


def main():
    parser = argparse.ArgumentParser(description=APP_NAME+" - "+APP_DESC+". Version " +
                                     APP_VERSION+" ("+APP_DATE+")",
                                     epilog="For more info visit "+APP_URL)
    parser.add_argument("-s", "--sample",
                        dest="sample_id",
                        default="",
                        help="Sample ID to phase (must exists in BAM and VCF)")
    parser.add_argument("-vi", "--vcf-in",
                        dest="vcf_in_path",
                        required=True,
                        help="Input VCF file")
    parser.add_argument("-vs", "--vcf-scaffold",
                        dest="vcf_scaffold_path",
                        default="",
                        help="Input phased files for scaffolding")
    parser.add_argument("-vo", "--vcf-out",
                        dest="vcf_out_path",
                        help="Output phased VCF file")
    parser.add_argument("-bi", "--bam-in",
                        dest="bam_in_path",
                        required=True,
                        help="Input BAM file")
    parser.add_argument("-bo", "--bam-out",
                        dest="bam_out_path",
                        default="",
                        help="Output BAM file (4 files will be created = " +
                        "one containing all read, 2 for each phase, " +
                        "and 1 for unphased reads)")
    parser.add_argument("-r", "--region",
                        dest="input_region",
                        default="",
                        help="Region specified as specified as: RNAME[:STARTPOS[-ENDPOS]]")
    parser.add_argument("-ei", "--bed-in",
                        dest="bed_in_path",
                        default="",
                        help="Input BED file with regions to phase. " +
                        "Used only if --region is not supplied")
    parser.add_argument("-eo", "--bed-out",
                        dest="bed_out_path",
                        default="",
                        help="Output BED file where haplotypes will be saved")
    parser.add_argument("-l", "--log-file",
                        dest="log_file_path",
                        default="",
                        help="Output log file")
    parser.add_argument("-no-pe", "--no-pair-end",
                        dest="pair_end",
                        action='store_false',
                        help="Don't use pair-end based haplotype merging")
    parser.add_argument("-fmq", "--filter-map-quality",
                        dest="filter_map_quality",
                        default=20,
                        type=int,
                        help="Minimum mapping quality to process")
    parser.add_argument("-md", "--max-depth",
                        dest="max_depth",
                        default=100,
                        type=int,
                        help="Maximum reads depth to examine")
    parser.add_argument("-q", "--quiet",
                        dest="quiet_mode",
                        action='store_true',
                        help="Quiet mode - don't output anything to console")

    args = vars(parser.parse_known_args()[0]).items()

    phaser = BAMPhaser()
    for param, value in args:
        setattr(phaser, param, value)

    phaser.phase()


if __name__ == "__main__":
    main()
