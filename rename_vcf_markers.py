import argparse
import gzip
import sys

def process_vcf(input_file, output_file, ref_index, ref_contigs, non_native_prefix=None):
    
    # Parse the comma-separated list of reference contigs into a set for fast lookup
    ref_contigs_set = set(ref_contigs.split(',')) if ref_contigs else set()
    
    # Determine the correct function to open the input file based on its extension
    opener = gzip.open if input_file.endswith('.gz') else open
    
    # Determine the correct function to open the output file based on its extension
    out_opener = gzip.open if output_file.endswith('.gz') else open
    
    # Open both the input file for reading (text mode) and output file for writing (text mode)
    with opener(input_file, 'rt') as f_in, out_opener(output_file, 'wt') as f_out:
        
        # Iterate through every line in the input file
        for line in f_in:
            
            # Check if the line is a header or metadata line
            if line.startswith('#'):
                
                # Write the header line to the output file directly
                f_out.write(line)
                
                # Move to the next line
                continue
            
            # Split the data line into its constituent columns using tab as delimiter
            parts = line.strip('\n').split('\t')
            
            # Extract the chromosome name from the first column
            chrom = parts[0]
            
            # Extract the physical position from the second column
            pos = parts[1]
            
            # Extract the current marker ID from the third column
            marker_id = parts[2]
            
            # Extract the INFO field from the eighth column
            info = parts[7]
            
            # Check if the marker ID is empty or missing (indicated by a period)
            if marker_id == '.':
                
                # Pad the physical position with leading zeros to ensure it is exactly 9 digits long
                pos_padded = str(pos).zfill(9)
                
                # Check if the chromosome belongs to the specified reference contigs
                if chrom in ref_contigs_set:
                    
                    # Strip 'Chr' from the chromosome name if it exists (e.g., 'Chr1A' becomes '1A')
                    chrom_suffix = chrom.replace('Chr', '') if chrom.startswith('Chr') else chrom
                    
                    # Construct the new marker ID using the reference index, stripped chromosome, and padded position
                    new_id = f"{ref_index}_{chrom_suffix}{pos_padded}"
                
                # Otherwise, if the chromosome does not belong to the reference contigs
                else:
                    
                    # Check if an optional non-native prefix was provided
                    if non_native_prefix:
                        
                        # Construct the ID using the contig name, the optional prefix, and padded position
                        new_id = f"{chrom}_{non_native_prefix}_{pos_padded}"
                        
                    # Otherwise, use the standard non-native format
                    else:
                        
                        # Construct the new marker ID using the full contig name and padded position
                        new_id = f"{chrom}_{pos_padded}"
                
                # Update the variant ID column with the newly constructed ID
                parts[2] = new_id
                
                # Split the INFO field by semicolons to check for the Original ID (OID) attribute
                info_parts = info.split(';')
                
                # Create a list to hold the updated INFO attributes
                new_info_parts = []
                
                # Iterate through each attribute in the INFO field
                for ip in info_parts:
                    
                    # Check if the attribute is an empty Original ID (OID=.)
                    if ip == 'OID=.':
                        
                        # Replace the empty OID with the new constructed ID
                        new_info_parts.append(f"OID={new_id}")
                        
                    # Otherwise, if it is not an empty OID
                    else:
                        
                        # Keep the attribute exactly as it was
                        new_info_parts.append(ip)
                
                # Rejoin the updated INFO attributes with semicolons and update the INFO column
                parts[7] = ';'.join(new_info_parts)
                
                # Join the updated columns back together with tabs to form the final line string
                line = '\t'.join(parts) + '\n'
            
            # Write the final line (whether updated or unchanged) to the output file
            f_out.write(line)

if __name__ == '__main__':
    
    # Initialize the argument parser for command line execution
    parser = argparse.ArgumentParser(description="Rename missing markers in VCF files.")
    
    # Add an argument for the input VCF file
    parser.add_argument('--input', required=True, help="Input VCF file (can be .gz)")
    
    # Add an argument for the output VCF file
    parser.add_argument('--output', required=True, help="Output VCF file (can be .gz)")
    
    # Add an argument for the reference index prefix
    parser.add_argument('--ref-index', required=True, help="Reference index string, e.g., CHS21")
    
    # Add an argument for the comma-separated list of reference contigs
    parser.add_argument('--ref-contigs', required=True, help="Comma separated list of contigs belonging to the reference")
    
    # Add an argument for an optional prefix on non-native contigs
    parser.add_argument('--non-native-prefix', required=False, default=None, help="Optional string to concatenate prior to position for non-native contigs")
    
    # Parse the arguments provided by the user
    args = parser.parse_args()
    
    # Run the main processing function with the provided arguments
    process_vcf(args.input, args.output, args.ref_index, args.ref_contigs, args.non_native_prefix)
