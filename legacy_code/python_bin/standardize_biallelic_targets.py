import argparse
import gzip
import os
import re
import sys
from collections import defaultdict


# Regular expression for valid canonical nucleotide sequences in VCF 4.2
VALID_BASES_RE = re.compile(r'^[ACGTNacgtn]+$')

# Regular expression for valid VCF 4.2 symbolic alleles (e.g., <INS>, <DEL> for presence/absence assays)
SYMBOLIC_ALLELE_RE = re.compile(r'^<[A-Za-z0-9_.:-]+>$')

# Regular expression to extract ID, Number, and Type from ##INFO and ##FORMAT VCF header lines
HEADER_META_RE = re.compile(r'^##(INFO|FORMAT)=<ID=([^,>]+),.*?Number=([^,>]+),Type=([^,>]+)')


def open_stream(file_path, mode):
    
    # Check if standard input or standard output was requested
    if file_path == '-':
        
        # Return standard input for reading or standard output for writing
        return sys.stdin if 'r' in mode else sys.stdout
    
    # Check if the file path ends with .gz for gzipped/BGZF stream handling
    if file_path.endswith('.gz'):
        
        # Open compressed stream in text mode
        return gzip.open(file_path, mode)
    
    # Open uncompressed plain-text file stream
    return open(file_path, mode, newline='\n' if 'w' in mode else None)


def parse_targets_manifest(targets_path, biallelic_mode):
    
    # Verify that the target manifest file exists on disk
    if not os.path.isfile(targets_path):
        sys.stderr.write(f"Error: Target manifest '{targets_path}' does not exist.\n")
        sys.exit(1)
    
    # Initialize list of parsed target dictionaries and lookup index by chromosome
    targets = []
    targets_by_chrom = defaultdict(list)
    seen_ids = set()
    
    # Open the manifest file and strip Windows/Unix line endings cleanly
    with open(targets_path, 'r', encoding='utf-8-sig') as f:
        raw_lines = [(idx + 1, line.rstrip('\r\n')) for idx, line in enumerate(f)]
    
    # Filter out completely empty or whitespace-only lines while retaining original line numbers
    non_empty = [(ln, txt) for ln, txt in raw_lines if txt.strip()]
    
    # Ensure the manifest contains at least one non-empty line
    if not non_empty:
        sys.stderr.write(f"Error: Target manifest '{targets_path}' is empty.\n")
        sys.exit(1)
    
    # Skip leading ## metadata comment lines if present
    start_idx = 0
    while start_idx < len(non_empty) and non_empty[start_idx][1].startswith('##'):
        start_idx += 1
    
    # Ensure rows remain after skipping ## comments
    if start_idx >= len(non_empty):
        sys.stderr.write(f"Error: Target manifest '{targets_path}' contains no header or target rows.\n")
        sys.exit(1)
    
    # Inspect the first non-metadata row to dynamically detect column names (tab-delimited or whitespace-delimited)
    first_ln, first_txt = non_empty[start_idx]
    first_cols = first_txt.split('\t') if '\t' in first_txt else first_txt.split()
    normalized_cols = [c.strip().lstrip('#').upper() for c in first_cols]
    
    # Determine whether the first row is a header row (starts with # or contains CHROM/POS/ID column names)
    has_header = (
        first_txt.startswith('#')
        or ('CHROM' in normalized_cols and 'POS' in normalized_cols)
        or ('SEQ_CHR' in normalized_cols and 'ID' in normalized_cols)
    )
    
    # Resolve column indices dynamically from header names if a header is present
    if has_header:
        col_map = {name: i for i, name in enumerate(normalized_cols)}
        
        # Locate CHROM, POS, and ID columns
        chrom_col = col_map.get('CHROM', col_map.get('CHR'))
        pos_col = col_map.get('POS')
        id_col = col_map.get('ID')
        
        # Verify required coordinate and ID columns exist
        missing_req = [name for name, idx in [('CHROM', chrom_col), ('POS', pos_col), ('ID', id_col)] if idx is None]
        if missing_req:
            sys.stderr.write(
                f"Error: Target manifest '{targets_path}' is missing required column(s): {', '.join(missing_req)}.\n"
                f"Detected header columns: {', '.join(first_cols)}\n"
            )
            sys.exit(1)
        
        # Locate REF and ALT columns when biallelic standardization mode is enabled
        ref_col = col_map.get('REF')
        alt_col = col_map.get('ALT')
        if biallelic_mode and (ref_col is None or alt_col is None):
            missing_ba = [name for name, idx in [('REF', ref_col), ('ALT', alt_col)] if idx is None]
            sys.stderr.write(
                f"Error: --biallelic-targets (-b) requires canonical REF and ALT columns in '{targets_path}', "
                f"but column(s) {', '.join(missing_ba)} were not found.\n"
                f"Detected header columns: {', '.join(first_cols)}\n"
            )
            sys.exit(1)
        
        # Data rows begin immediately after the header row
        data_rows = non_empty[start_idx + 1:]
    
    # Otherwise fall back to standard positional columns (#CHROM=0, POS=1, ID=2, [REF=3, ALT=4])
    else:
        if len(first_cols) < (5 if biallelic_mode else 3):
            req_n = 5 if biallelic_mode else 3
            sys.stderr.write(
                f"Error: Headerless target manifest '{targets_path}' has {len(first_cols)} column(s) at line {first_ln}, "
                f"but at least {req_n} tab-delimited columns are required.\n"
            )
            sys.exit(1)
        chrom_col, pos_col, id_col = 0, 1, 2
        ref_col, alt_col = (3, 4) if len(first_cols) >= 5 else (None, None)
        data_rows = non_empty[start_idx:]
    
    # Ensure at least one target locus row is present
    if not data_rows:
        sys.stderr.write(f"Error: Target manifest '{targets_path}' contains a header but zero target rows.\n")
        sys.exit(1)
    
    # Determine minimum required column index per row
    used_cols = [chrom_col, pos_col, id_col] + ([ref_col, alt_col] if biallelic_mode else [])
    min_cols_needed = max(used_cols) + 1
    
    # Parse and validate every target row
    for ln, txt in data_rows:
        
        # Skip inline comment lines if any exist
        if txt.startswith('#'):
            continue
        
        cols = txt.split('\t') if '\t' in txt else txt.split()
        if len(cols) < min_cols_needed:
            sys.stderr.write(
                f"Error: Line {ln} in '{targets_path}' has {len(cols)} tab-delimited column(s), "
                f"but at least {min_cols_needed} are required.\n"
            )
            sys.exit(1)
        
        chrom = cols[chrom_col].strip()
        pos_str = cols[pos_col].strip()
        tid = cols[id_col].strip()
        
        # Validate non-empty CHROM and ID
        if not chrom or not tid:
            sys.stderr.write(f"Error: Empty CHROM or ID value at line {ln} in '{targets_path}'.\n")
            sys.exit(1)
        
        # Validate positive integer POS
        if not pos_str.isdigit() or int(pos_str) <= 0:
            sys.stderr.write(
                f"Error: Invalid genomic coordinate '{pos_str}' at line {ln} in '{targets_path}'. "
                f"POS must be a positive integer.\n"
            )
            sys.exit(1)
        pos = int(pos_str)
        
        # Parse and validate REF and ALT when biallelic mode is active
        ref_val = None
        alt_val = None
        symbolic_assay = False
        if biallelic_mode:
            ref_val = cols[ref_col].strip().upper()
            alt_val = cols[alt_col].strip()
            
            # Normalize standard nucleotide ALT to uppercase while preserving case of symbolic tags like <INS>/<DEL>
            if not SYMBOLIC_ALLELE_RE.match(alt_val):
                alt_val = alt_val.upper()
            else:
                symbolic_assay = True
            
            # Disallow multiallelic comma-separated target ALTs in strict biallelic manifest
            if ',' in alt_val or ',' in ref_val:
                sys.stderr.write(
                    f"Error: Multiallelic manifest entry at line {ln} ({tid}: REF='{ref_val}', ALT='{alt_val}'). "
                    f"Each target must specify a single canonical REF and ALT allele.\n"
                )
                sys.exit(1)
            
            # Validate REF nucleotide characters
            if not VALID_BASES_RE.match(ref_val):
                sys.stderr.write(
                    f"Error: Invalid REF allele '{ref_val}' for marker '{tid}' at line {ln} in '{targets_path}'.\n"
                )
                sys.exit(1)
            
            # Validate ALT nucleotide characters or VCF 4.2 symbolic allele
            if not (VALID_BASES_RE.match(alt_val) or symbolic_assay):
                sys.stderr.write(
                    f"Error: Invalid ALT allele '{alt_val}' for marker '{tid}' at line {ln} in '{targets_path}'.\n"
                )
                sys.exit(1)
        
        # Skip exact duplicate marker IDs while keeping first canonical definition
        if tid in seen_ids:
            continue
        seen_ids.add(tid)
        
        target_entry = {
            'order': len(targets),
            'chrom': chrom,
            'pos': pos,
            'id': tid,
            'ref': ref_val,
            'alt': alt_val,
            'symbolic': symbolic_assay,
        }
        targets.append(target_entry)
        targets_by_chrom[chrom].append(target_entry)
    
    return targets, targets_by_chrom


def extract_id_and_oid_set(vid_col, info_col):
    
    # Collect all semicolon-delimited IDs from VCF column 3
    ids = set(x for x in vid_col.split(';') if x and x != '.')
    
    # Also collect all comma-delimited Hotspot OIDs from INFO/OID (survives AgriSeq_Merger.sh -d)
    for token in info_col.split(';'):
        if token.startswith('OID='):
            for oid in token[4:].split(','):
                if oid and oid != '.':
                    ids.add(oid)
            break
    
    return ids


def evaluate_candidate_match(target, pos_vcf, ref_vcf, orig_alts, biallelic_mode):
    
    # Extract target coordinate
    t_pos = target['pos']
    
    # In non-biallelic default subsetting mode, match exact POS or 5'-spanning haplotype block
    if not biallelic_mode:
        if pos_vcf == t_pos:
            return (1, 0, None)
        if pos_vcf < t_pos and (pos_vcf + len(ref_vcf) - 1) >= t_pos:
            return (2, t_pos - pos_vcf, None)
        return None
    
    # Symbolic presence/absence assays (<INS>, <DEL>) cannot match nucleotide variant calls
    if target['symbolic']:
        return None
    
    t_ref = target['ref']
    t_alt = target['alt']
    ref_u = ref_vcf.upper()
    alts_u = [a.upper() for a in orig_alts]
    
    # Tier 1: Exact start coordinate match where VCF REF equals canonical REF + 3' padding suffix
    if pos_vcf == t_pos and ref_u.startswith(t_ref):
        suffix = ref_u[len(t_ref):]
        expected_alt = t_alt + suffix
        alt_idx = alts_u.index(expected_alt) if expected_alt in alts_u else None
        return (1, 0, alt_idx)
    
    # Tier 2: 5'-anchored TVC Hotspot block spanning across canonical POS (e.g., 1-bp MNP anchor trim)
    if pos_vcf < t_pos and (pos_vcf + len(ref_u) - 1) >= t_pos:
        offset = t_pos - pos_vcf
        prefix = ref_u[:offset]
        rem_ref = ref_u[offset:]
        if rem_ref.startswith(t_ref):
            suffix = rem_ref[len(t_ref):]
            expected_alt = prefix + t_alt + suffix
            alt_idx = alts_u.index(expected_alt) if expected_alt in alts_u else None
            return (2, offset, alt_idx)
    
    # Tier 3: Equivalent 1-bp right-shifted homopolymer/repeat indel fallback (only used if Tier 1/2 absent)
    if pos_vcf == t_pos + 1:
        matched_idx = None
        for idx, cand_alt in enumerate(alts_u):
            if t_ref + cand_alt == t_alt + ref_u:
                matched_idx = idx
                break
        if matched_idx is not None:
            return (3, -1, matched_idx)
    
    return None


def slice_meta_value(val, num_spec, val_type, n_alts, target_alt_idx, is_called=True, tag_id=''):
    
    # Preserve missing scalar values directly
    if val == '.' or val == '':
        return '.'
    
    # Handle Number=A (one value per ALT allele)
    if num_spec == 'A':
        parts = val.split(',')
        if target_alt_idx is not None and target_alt_idx < len(parts):
            return parts[target_alt_idx]
        if val_type in ('Integer', 'Float') and is_called:
            return '0'
        return '.'
    
    # Handle Number=R (one value per allele including REF, i.e., REF + each ALT)
    if num_spec == 'R':
        parts = val.split(',')
        ref_part = parts[0] if parts else '.'
        if target_alt_idx is not None and (1 + target_alt_idx) < len(parts):
            alt_part = parts[1 + target_alt_idx]
        elif val_type in ('Integer', 'Float') and is_called and ref_part != '.':
            alt_part = '0'
        else:
            alt_part = '.'
        return f"{ref_part},{alt_part}"
    
    # Handle Number=G (one value per diploid genotype: 0/0, 0/1, 1/1, ...)
    if num_spec == 'G':
        parts = val.split(',')
        if target_alt_idx is not None:
            t = 1 + target_alt_idx
            i00 = 0
            i01 = t * (t + 1) // 2
            i11 = t * (t + 1) // 2 + t
            if i11 < len(parts):
                return f"{parts[i00]},{parts[i01]},{parts[i11]}"
        return '.'
    
    # Handle TVC's INFO/FR tag (declared Number=. in TVC header but populated per ALT allele)
    if tag_id == 'FR' and n_alts > 1 and ',' in val:
        parts = val.split(',')
        if len(parts) == n_alts:
            if target_alt_idx is not None and target_alt_idx < len(parts):
                return parts[target_alt_idx]
            return '.'
    
    # Return all other tags (Number=0, Number=1, etc.) unchanged
    return val


def infer_canonical_type_and_len(ref_seq, alt_seq):
    
    # Determine canonical TVC variant TYPE and LEN for the standardized biallelic REF and ALT
    r_len = len(ref_seq)
    a_len = len(alt_seq)
    if r_len == 1 and a_len == 1:
        return 'snp', 1
    if r_len == a_len:
        return 'mnp', a_len
    if a_len > r_len and (alt_seq.startswith(ref_seq) or alt_seq.endswith(ref_seq)):
        return 'ins', a_len - r_len
    if r_len > a_len and (ref_seq.startswith(alt_seq) or ref_seq.endswith(alt_seq)):
        return 'del', r_len - a_len
    return 'complex', a_len


def standardize_info_field(info_str, info_meta, n_alts, target_alt_idx, target, stats=None, has_samples=False):
    
    # Preserve empty/missing INFO column as '.'
    if not info_str or info_str == '.':
        return '.'
    
    # Pre-compute canonical variant TYPE and LEN for the target biallelic alleles
    canon_type, canon_len = infer_canonical_type_and_len(target['ref'], target['alt'])
    
    out_tokens = []
    for token in info_str.split(';'):
        if not token:
            continue
        if '=' not in token:
            out_tokens.append(token)
            continue
        
        k, v = token.split('=', 1)
        
        # Canonicalize TVC Hotspot provenance tags to match the standardized target locus
        if k == 'OID':
            out_tokens.append(f"OID={target['id']}")
            continue
        if k == 'OPOS':
            out_tokens.append(f"OPOS={target['pos']}")
            continue
        if k == 'OREF':
            out_tokens.append(f"OREF={target['ref']}")
            continue
        if k == 'OALT':
            out_tokens.append(f"OALT={target['alt']}")
            continue
        if k == 'OMAPALT':
            out_tokens.append(f"OMAPALT={target['alt']}")
            continue
        if k == 'TYPE':
            out_tokens.append(f"TYPE={canon_type}")
            continue
        if k == 'LEN':
            out_tokens.append(f"LEN={canon_len}")
            continue
        
        # Update cohort genotype and allele count summary tags when sample columns are present
        if has_samples and stats is not None:
            if k == 'AC':
                out_tokens.append(f"AC={stats['ac']}")
                continue
            if k == 'AN':
                out_tokens.append(f"AN={stats['an']}")
                continue
            if k == 'NS':
                out_tokens.append(f"NS={stats['called']}")
                continue
            if k == 'WT':
                out_tokens.append(f"WT={stats['hom_ref']}")
                continue
            if k == 'HET':
                out_tokens.append(f"HET={stats['het']}")
                continue
            if k == 'HOM':
                out_tokens.append(f"HOM={stats['hom_alt']}")
                continue
            if k == 'NC':
                out_tokens.append(f"NC={stats['missing']}")
                continue
        
        num_spec, val_type = info_meta.get(k, ('1', 'String'))
        sliced_v = slice_meta_value(v, num_spec, val_type, n_alts, target_alt_idx, is_called=True, tag_id=k)
        out_tokens.append(f"{k}={sliced_v}")
    
    return ';'.join(out_tokens) if out_tokens else '.'


def standardize_record(parts, target, target_alt_idx, info_meta, format_meta):
    
    # Record original VCF coordinates and alleles for QC reporting
    orig_pos = int(parts[1])
    orig_ref = parts[3]
    orig_alt_str = parts[4]
    orig_alts = [] if orig_alt_str == '.' else orig_alt_str.split(',')
    n_alts = len(orig_alts)
    
    # Determine whether coordinates or alleles required standardization
    was_modified = (
        orig_pos != target['pos']
        or parts[2] != target['id']
        or orig_ref.upper() != target['ref']
        or orig_alt_str.upper() != target['alt']
    )
    
    # Rewrite canonical VCF locus columns (#CHROM, POS, ID, REF, ALT)
    parts[0] = target['chrom']
    parts[1] = str(target['pos'])
    parts[2] = target['id']
    parts[3] = target['ref']
    parts[4] = target['alt']
    
    # Initialize genotype and allele counters for QC audit reporting and INFO summary updates
    stats = {
        'called': 0,
        'hom_ref': 0,
        'het': 0,
        'hom_alt': 0,
        'masked_offtarget': 0,
        'missing': 0,
        'ac': 0,
        'an': 0,
    }
    has_samples = len(parts) > 8
    
    # If VCF has FORMAT and sample columns, remap genotypes and slice Number=A/R/G FORMAT fields
    if has_samples:
        fmt_keys = parts[8].split(':')
        fmt_specs = [format_meta.get(k, ('1', 'String')) for k in fmt_keys]
        
        # Build allele index lookup mapping strictly '0' -> '0' and target ALT -> '1'
        valid_allele_map = {'0': '0'}
        if target_alt_idx is not None:
            valid_allele_map[str(1 + target_alt_idx)] = '1'
        
        # Pre-compute clean missing sample template (./.:.:...) matching FORMAT column count
        missing_tail = [':'.join(['.'] * (len(fmt_keys) - 1))] if len(fmt_keys) > 1 else []
        
        for s_idx in range(9, len(parts)):
            s_val = parts[s_idx]
            subfields = s_val.split(':')
            if len(subfields) < len(fmt_keys):
                subfields.extend(['.'] * (len(fmt_keys) - len(subfields)))
            
            gt_raw = subfields[0]
            sep = '|' if '|' in gt_raw else ('/' if '/' in gt_raw else None)
            gt_tokens = re.split(r'[|/]', gt_raw) if sep else [gt_raw]
            
            # Case 1: Sample was already completely missing (./. or .)
            if all(tok == '.' for tok in gt_tokens):
                stats['missing'] += 1
                new_sub = [gt_raw]
                for f_i in range(1, len(fmt_keys)):
                    num_spec, val_type = fmt_specs[f_i]
                    orig_sub = subfields[f_i]
                    # Preserve 0s in TVC's ./.:0:0:... missing records when slicing Number=A/R
                    treat_zero = (orig_sub == '0' or orig_sub.startswith('0,'))
                    new_sub.append(
                        slice_meta_value(orig_sub, num_spec, val_type, n_alts, target_alt_idx, is_called=treat_zero, tag_id=fmt_keys[f_i])
                    )
                parts[s_idx] = ':'.join(new_sub)
                continue
            
            # Case 2: Sample carries an off-target allele or half-missing call -> mask completely to ./.
            if any(tok not in valid_allele_map for tok in gt_tokens):
                stats['masked_offtarget'] += 1
                stats['missing'] += 1
                masked_gt = './.' if len(gt_tokens) >= 2 else '.'
                parts[s_idx] = ':'.join([masked_gt] + missing_tail)
                was_modified = True
                continue
            
            # Case 3: Valid biallelic target genotype (composed exclusively of canonical REF and/or ALT)
            mapped = [valid_allele_map[tok] for tok in gt_tokens]
            
            # Canonicalize unphased heterozygous allele order (1/0 -> 0/1)
            if sep == '/' and len(mapped) == 2 and mapped[0] == '1' and mapped[1] == '0':
                mapped = ['0', '1']
            
            new_gt = sep.join(mapped) if sep else mapped[0]
            if new_gt != gt_raw:
                was_modified = True
            
            # Tally genotype and allele counts
            stats['called'] += 1
            stats['an'] += len(mapped)
            stats['ac'] += sum(1 for m in mapped if m == '1')
            if all(m == '0' for m in mapped):
                stats['hom_ref'] += 1
            elif all(m == '1' for m in mapped):
                stats['hom_alt'] += 1
            else:
                stats['het'] += 1
            
            # Slice per-allele FORMAT fields to the target ALT allele
            new_sub = [new_gt]
            for f_i in range(1, len(fmt_keys)):
                num_spec, val_type = fmt_specs[f_i]
                new_sub.append(
                    slice_meta_value(subfields[f_i], num_spec, val_type, n_alts, target_alt_idx, is_called=True, tag_id=fmt_keys[f_i])
                )
            parts[s_idx] = ':'.join(new_sub)
    
    # Slice per-allele INFO attributes and update cohort counts cleanly
    parts[7] = standardize_info_field(
        parts[7], info_meta, n_alts, target_alt_idx, target, stats=stats, has_samples=has_samples
    )
    
    return parts, orig_pos, orig_ref, orig_alt_str, was_modified, stats


def process_vcf(input_vcf, output_vcf, targets_path, biallelic_mode=False, qc_output_path=None):
    
    # Parse and validate the target manifest dynamically
    targets, targets_by_chrom = parse_targets_manifest(targets_path, biallelic_mode)
    
    # Storage for VCF header metadata and best matched record per target ID
    header_lines = []
    info_meta = {}
    format_meta = {}
    best_match_by_id = {}
    
    # Stream through the input VCF to parse headers and match target records
    with open_stream(input_vcf, 'rt') as f_in:
        for line in f_in:
            if line.startswith('#'):
                header_lines.append(line)
                m = HEADER_META_RE.match(line)
                if m:
                    sec, tag_id, num_spec, val_type = m.groups()
                    if sec == 'INFO':
                        info_meta[tag_id] = (num_spec, val_type)
                    else:
                        format_meta[tag_id] = (num_spec, val_type)
                continue
            
            parts = line.rstrip('\r\n').split('\t')
            if len(parts) < 8:
                continue
            
            chrom = parts[0]
            if chrom not in targets_by_chrom:
                continue
            
            try:
                pos_vcf = int(parts[1])
            except ValueError:
                continue
            
            vid_col = parts[2]
            ref_vcf = parts[3]
            alt_col = parts[4]
            info_col = parts[7]
            orig_alts = [] if alt_col == '.' else alt_col.split(',')
            
            # Extract all marker IDs associated with this VCF row (column 3 ID + INFO/OID)
            row_ids = extract_id_and_oid_set(vid_col, info_col)
            if not row_ids:
                continue
            
            # Test candidate targets on this chromosome whose ID is present on this row
            for target in targets_by_chrom[chrom]:
                tid = target['id']
                if tid not in row_ids:
                    continue
                
                res = evaluate_candidate_match(target, pos_vcf, ref_vcf, orig_alts, biallelic_mode)
                if res is None:
                    continue
                
                tier, offset, alt_idx = res
                
                # Rank candidate matches: lower tier wins (Tier 1 exact > Tier 2 5'-trim > Tier 3 1bp-shift),
                # breaking ties by presence of target ALT in orig_alts
                rank_key = (tier, 0 if alt_idx is not None else 1, abs(offset))
                prev = best_match_by_id.get(tid)
                if prev is None or rank_key < prev['rank']:
                    best_match_by_id[tid] = {
                        'rank': rank_key,
                        'tier': tier,
                        'alt_idx': alt_idx,
                        'parts': list(parts),
                    }
    
    # Open output VCF stream and write header + matched/standardized records in manifest order
    qc_rows = []
    matched_loci = 0
    modified_loci = 0
    total_masked_genotypes = 0
    emitted_raw_keys = set()
    
    with open_stream(output_vcf, 'wt') as f_out:
        for hline in header_lines:
            f_out.write(hline if hline.endswith('\n') else hline + '\n')
        
        for target in targets:
            tid = target['id']
            match = best_match_by_id.get(tid)
            
            # Handle targets absent from the input VCF
            if match is None:
                if biallelic_mode:
                    status = 'SYMBOLIC_UNCALLED_ASSAY' if target['symbolic'] else 'MISSING_FROM_VCF'
                    qc_rows.append((
                        target['chrom'], target['pos'], tid, target['ref'], target['alt'],
                        status, '.', '.', '.', '.', 0, 0, 0, 0, 0, 0
                    ))
                continue
            
            matched_loci += 1
            
            # Default non-destructive mode: write matched VCF record completely untouched (deduplicating shared blocks)
            if not biallelic_mode:
                raw_key = tuple(match['parts'][:5])
                if raw_key not in emitted_raw_keys:
                    emitted_raw_keys.add(raw_key)
                    f_out.write('\t'.join(match['parts']) + '\n')
                continue
            
            # Strict biallelic mode: rewrite coordinates/alleles, slice INFO/FORMAT, and mask off-target GTs
            std_parts, orig_pos, orig_ref, orig_alt_str, was_mod, stats = standardize_record(
                list(match['parts']), target, match['alt_idx'], info_meta, format_meta
            )
            f_out.write('\t'.join(std_parts) + '\n')
            
            if was_mod:
                modified_loci += 1
            total_masked_genotypes += stats['masked_offtarget']
            
            status = 'STANDARDIZED' if was_mod else 'CLEAN_BIALLELIC'
            alt_idx_str = str(match['alt_idx'] + 1) if match['alt_idx'] is not None else 'NONE'
            qc_rows.append((
                target['chrom'], target['pos'], tid, target['ref'], target['alt'],
                status, str(orig_pos), orig_ref, orig_alt_str, alt_idx_str,
                stats['called'], stats['hom_ref'], stats['het'], stats['hom_alt'],
                stats['masked_offtarget'], stats['missing']
            ))
    
    # Write companion per-marker QC audit TSV if requested in biallelic mode
    if biallelic_mode and qc_output_path:
        qc_dir = os.path.dirname(os.path.abspath(qc_output_path))
        if qc_dir and not os.path.isdir(qc_dir):
            os.makedirs(qc_dir, exist_ok=True)
        with open(qc_output_path, 'w', newline='\n') as qf:
            qf.write(
                "#CHROM\tPOS\tID\tREF\tALT\tSTATUS\tORIGINAL_POS\tORIGINAL_REF\tORIGINAL_ALT\t"
                "TARGET_ALT_INDEX\tCALLED_GENOTYPES\tHOM_REF_COUNT\tHET_COUNT\tHOM_ALT_COUNT\t"
                "MASKED_OFFTARGET_GENOTYPES\tMISSING_GENOTYPES\n"
            )
            for row in qc_rows:
                qf.write('\t'.join(str(x) for x in row) + '\n')
    
    # Print concise execution summary to stderr/stdout (stderr when output_vcf is '-', else stdout)
    log_stream = sys.stderr if output_vcf == '-' else sys.stdout
    missing_loci = len(targets) - matched_loci
    if biallelic_mode:
        log_stream.write(f"  -> Biallelic Target Standardization Summary:\n")
        log_stream.write(f"       Manifest targets evaluated : {len(targets)}\n")
        log_stream.write(f"       Target loci matched in VCF : {matched_loci}\n")
        log_stream.write(f"       Loci standardized/trimmed  : {modified_loci}\n")
        log_stream.write(f"       Off-target GTs masked (./.): {total_masked_genotypes}\n")
        log_stream.write(f"       Targets absent from VCF    : {missing_loci}\n")
        if qc_output_path:
            log_stream.write(f"       Companion QC report        : {qc_output_path}\n")
    else:
        log_stream.write(
            f"  -> Coordinate + ID verification retained {matched_loci} of {len(targets)} manifest targets "
            f"({missing_loci} absent from VCF).\n"
        )


if __name__ == '__main__':
    
    # Initialize command-line argument parser
    parser = argparse.ArgumentParser(
        description=(
            "Filter VCF records by manifest coordinate + ID, and optionally enforce strict "
            "canonical biallelic REF/ALT alleles (masking off-target calls to ./.) for T3/GWAS/GS."
        )
    )
    
    # Add argument for input VCF file
    parser.add_argument('--input', required=True, help="Input VCF path (.vcf, .vcf.gz, or '-' for stdin)")
    
    # Add argument for output VCF file
    parser.add_argument('--output', required=True, help="Output VCF path (.vcf, .vcf.gz, or '-' for stdout)")
    
    # Add argument for target manifest TSV
    parser.add_argument(
        '--targets',
        required=True,
        help="Tab-delimited target manifest (#CHROM POS ID [REF ALT])"
    )
    
    # Add optional flag to activate strict biallelic standardization
    parser.add_argument(
        '-b', '--biallelic',
        action='store_true',
        help="Enforce canonical REF/ALT from manifest, trim padded windows, slice Number=A/R/G tags, and mask off-target genotypes to ./."
    )
    
    # Add optional companion QC TSV output path
    parser.add_argument(
        '--qc-output',
        required=False,
        default=None,
        help="Optional output path for per-marker biallelic QC summary TSV"
    )
    
    # Parse CLI arguments and run VCF processing
    args = parser.parse_args()
    process_vcf(
        input_vcf=args.input,
        output_vcf=args.output,
        targets_path=args.targets,
        biallelic_mode=args.biallelic,
        qc_output_path=args.qc_output,
    )
