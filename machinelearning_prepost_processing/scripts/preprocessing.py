"""preprocessing.py: Data preprocessing
Adapated from CNN pipeline mouse_sst/preprocessing.py

Usage:
python preprocessing.py expand_peaks -i <input bed files> -o <output bed files> -l 500
"""
import operator
import os

import numpy as np
import pandas as pd
from pybedtools import BedTool
from tqdm import tqdm


SPLIT_REFERENCE_SPECIES_MAPPING = {
	'human': 'hg38',
	'homo_sapiens': 'hg38',
	'hg38': 'hg38',
	'mouse': 'mm10',
	'mus_musculus': 'mm10',
	'mm10': 'mm10',
	'rhesus': 'rheMac8',
	'macaca_mulatta': 'rheMac8',
	'rhemac8': 'rheMac8',
	'rhemac10': 'rheMac10',
	'rat': 'rn6',
	'rattus_norvegicus': 'rn6',
	'rn6': 'rn6'
}


def main(args):
	if args.function == 'expand_peaks':
		expand_peaks(args.in_files, args.out_files, args.length)
	elif args.function == 'subtract':
		subtract(args.include, args.exclude, out_file=args.out_files, window=args.window)
	elif args.function == 'get_original':
		get_original(args.in_files, args.out_files)
	elif args.function == 'filter_length':
		filter_length(args.in_files, args.out_files, args.conditions, combination=args.combination)
	elif args.function == 'filter_chromosome':
		filter_chromosome(args.in_files, args.out_files, drop_non_chromosomal=args.drop_non_chromosomal,
			drop_values=args.drop_values)
	elif args.function == 'split':
		split(args.in_files, args.split_ref_species, args.out_dir,
			split_ref_dir=args.split_ref_dir, do_get_original=args.get_original,
			create_folds=args.create_folds)
	else:
		raise ValueError(f"Invalid args: {args}")

def expand_peaks(in_files, out_files, length):
	"""Standardize bed, narrowPeak, or summits.rds file peaks to a uniform length.
	Adapted from expand_peaks.py by Calvin Chen.

	Actions:
	- Expand peaks to the same length
	- Drop duplicate peaks

	Args:
		in_files (list of str): Can be .narrowPeak, .narrowPeak.gz, or summits.rds
		out_files (list of str): Can be .narrowPeak or .narrowPeak.gz
		length (int): Length to expand peaks to.
	"""
	# Check args
	_check_parallel_files(in_files, out_files)

	if length is None:
		raise ValueError(f"Missing required parameter `length` in expand_peaks")

	# Expand peaks
	messages = []
	for in_file, out_file in tqdm(zip(in_files, out_files), total=len(in_files)):
		# Skip if file exists
		if os.path.exists(out_file):
			print(f"Output file {out_file} exists, skipping")
			continue

		# Read input file
		if in_file.endswith('.rds'):
			# Convert GenomicRanges object to df
			# NOTE: It would be faster to call read_gr_rds() only once, with a list of files,
			# because of the upfront cost of loading the GenomicRanges library.
			# Deciding to leave it like this for now because it's simpler.
			df = read_gr_rds(in_file)
			origin = 'rds'
		else:
			df = pd.read_csv(in_file, delim_whitespace=True, header=None)
			origin = 'narrowPeak'
		orig_len = len(df)

		# Compute peak boundaries and summits
		df = _expand(df, length, origin)

		# Drop duplicate peaks
		df = _drop_bed_duplicates(df)

		# Make directory and write result .bed or .narrowPeak
		_make_parent_dir(out_file)
		df.to_csv(out_file, index=False, sep="\t", header=None)

		# Update report
		messages.append(f"{out_file}: Expanded {len(df)} of {orig_len} peaks.")

	# Report results
	for message in messages:
		print(message)

def _expand(df, length, origin):
	"""Compute peak boundaries and summits."""
	if origin == 'narrowPeak':
		# Assumes this df came from a narrowPeak file. Columns (0-based):
		# 1: left endpoint.
		# 2: right endpoint.
		# 9: summit offset.

		# Compute summit index
		summits = df[1] + df[9]
		# Compute boundaries
		df = _compute_endpoints(df, summits, length)
		# Recompute summit offset
		df[9] = summits - df[1]

	elif origin == 'rds':
		# Assumes this df came from a `summits.rds` file. Bed format with columns (0-based):
		# 1: left endpoint. Since this is from summits.rds, we assume that this is the summit.
		# 2: right endpoint.
		# 4: signal
		
		# Get summit as left endpoint
		summits = df[1]
		# Compute boundaries
		df = _compute_endpoints(df, summits, length)

	else:
		raise NotImplementedError("Invalid origin")

	return df

def _compute_endpoints(df, summits, length):
	"""Compute left and right endpoints by expanding equally from summit"""
	df[1] = np.floor(summits - length / 2).astype(int)
	df[2] = df[1] + length
	return df

def read_gr_rds(path):
	"""Read GenomicRanges object from .rds file and convert to pd DataFrame.

	NOTE: This implementation is kind of silly because we read a file, write another, and read again.
	Technically we should be able to work with the GenomicRanges object directly using rpy2.
	This alternative solution was chosen because it's less complex.
	"""
	import tempfile
	import subprocess
	with tempfile.NamedTemporaryFile(suffix='.bed') as fp:
		subprocess.call(["Rscript", "--vanilla", "gr_to_bed.r", path, fp.name])
		bed = pd.read_csv(fp.name, delim_whitespace=True, header=None)
	return bed

def subtract(include_files, exclude_files, out_files=None, window=0):
	"""Take the union of included peaks, and subtract excluded peaks.

	Args:
		include_files (list of str): paths to .bed or .narrowPeak files of
			peaks to include.
		exclude_files (list of str): paths to .bed or .narrowPeak files of
			peaks to exclude. The peaks from include_files will only be included if
			they do not overlap any peaks from any of these files.
		out_files (str, optional): path to .bed or .narrowPeak file to save the result.
		window (int, default 0): the same as the -w flag for bedtools window.
			If window > 0, then this expands the region that is considered overlapping.
	"""
	# Handle arguments
	if out_files is not None:
		if isinstance(out_files, list):
			if len(out_files) > 1:
				raise ValueError("Multiple output files given, please pass only one output file")
			out_files = out_files[0]
		if not isinstance(out_files, str):
			raise ValueError(f"Invalid type for out_file, expected str but got {type(out_files)}, {out_files}")
	out_file = out_files

	# Union
	dfs = [pd.read_csv(in_file, delim_whitespace=True, header=None, na_values=['.'])
			for in_file in include_files]
	df = pd.concat(dfs)
	orig_len = len(df)
	# Drop duplicates
	df = _drop_bed_duplicates(df)
	# Replace NaN values with "." for a well-formed .bed file
	df = df.replace(np.NaN, '.')
	cur_len = len(df)
	if cur_len < orig_len:
		print(f"Info (union): Dropped {orig_len - cur_len} duplicates; kept {cur_len} of {orig_len}")

	# Convert df to BedTool
	bt = BedTool([[str(itm) for itm in row] for row in df.values.tolist()])

	# Subtraction
	if exclude_files:
		orig_len = len(bt)
		for file in exclude_files:
			bt = bt.window(BedTool(file), v=True, w=window)
		cur_len = len(bt)
		print(f"Info (subtract): Subtracted {orig_len - cur_len} peaks; kept {cur_len} of {orig_len}")

	# Save to file
	if out_file is not None:
		_make_parent_dir(out_file)
		bt.saveas(out_file)

	return bt

def get_original(in_files, out_files):
	"""Get original coordinates of mapped peaks from the name column.

	Args:
		in_files (list of str): can be .bed, .narrowPeak, .bed.gz, or .narrowPeak.gz
			name column must be [genome (optional)]:[chr]:[start]-[end]
		out_files (list of str): can be .bed, .narrowPeak, .bed.gz, or .narrowPeak.gz
	"""
	# Check args
	_check_parallel_files(in_files, out_files)

	# Get original coordinates
	messages = []
	for in_file, out_file in zip(in_files, out_files):
		df = pd.read_csv(in_file, delim_whitespace=True, header=None)
		orig_len = len(df)

		# Get chr_name, start, and end from name column 3 (0-based)
		ranges = pd.DataFrame(
			df[3].apply(_get_range).tolist(),
			columns=['genome', 'chr', 'start', 'end', 'summit']
			# Keep int columns as int, if they had any N/A values
			).convert_dtypes()

		# Re-assign chr_name, start, and end
		df[0] = ranges['chr']
		df[1] = ranges['start']
		df[2] = ranges['end']
		# Re-assign summit offset if it exists and if the name column had at least one not-null summit value
		if (9 in df.columns) and (ranges['summit'].notnull().any()):
			# Replace any null values with -1
			ranges['summit'] = ranges['summit'].fillna(-1)
			df[9] = ranges['summit']

		# Drop rows if the new chr_name, start, or end are null
		df = df.dropna(subset=[0, 1, 2])

		messages.append(f"Info (get_original): Recovered {len(df)} of {orig_len} peaks from {in_file}")

		# Make directory and write result .bed or .narrowPeak
		_make_parent_dir(out_file)
		df.to_csv(out_file, index=False, sep="\t", header=None)

	# Report
	for message in messages:
		print(message)

def _get_range_old(name):
	"""From [genome (optional)]:[chr]:[start]-[end], extract genome, chr, start, and end

	Args:
		name (str): expected to be [genome (optional)]:[chr]:[start]-[end]

	Returns: tuple
		genome (str)
		chr_name (str)
		start (int)
		end (int)
	"""
	# Get name parts as list split by ':'
	name = name.split(':')
	if (len(name) < 2) or (len(name) > 3):
		return None, None, None, None

	# Get start and end
	start_end = name[-1].split('-')
	if len(start_end) != 2:
		return None, None, None, None
	try:
		start = int(start_end[0])
		end = int(start_end[1])
	except ValueError:
		return None, None, None, None

	chr_name = name[-2]
	genome = None
	if len(name) > 2:
		genome = name[-3]

	return genome, chr_name, start, end

def _get_range(name):
	"""From [genome(optional)]:[chr]:[start]-[end]:[summit (optional)],
	extract genome, chr, start, end, and summit

	Args:
		name (str): expected to be [genome (optional)]:[chr]:[start]-[end]:[summit (optional)]

	Returns: tuple
		genome (str)
		chr_name (str)
		start (int)
		end (int)
		summit (int)
	"""
	# Get name parts as list split by ':'
	parts = name.split(':')
	if (len(parts) < 2) or (len(parts) > 4):
		print(f"Warning (_get_range): invalid ':'. Got '{name}'")
		return None, None, None, None, None

	# Find the unique part that contains '-' and assume this is [start]-[end]
	idxes = [i for i in range(len(parts)) if '-' in parts[i]]
	if len(idxes) != 1:
		print(f"Warning (_get_range): cannot parse [start]-[end]. Got '{name}'")
		return None, None, None, None, None
	start_end_idx = idxes[0]
	# Check that this part is either second or third
	if start_end_idx not in [1, 2]:
		print(f"Warning (_get_range): [start]-[end] not in expected place. Got '{name}'")
		return None, None, None, None, None
	# Check that this part contains exactly one '-' and split into start and end
	start_end_parts = parts[start_end_idx].split('-')
	if len(start_end_parts) != 2:
		print(f"Warning (_get_range): [start]-[end] parsed in an unexpected way. Got '{name}'")
		return None, None, None, None, None
	start, end = start_end_parts[0], start_end_parts[1]
	# Check that start and end are int
	try:
		start, end = int(start), int(end)
	except Exception as e:
		print(f"Warning (_get_range): start and end are not integers. Got '{name}'")
		return None, None, None, None, None

	# Get genome
	genome = None
	if start_end_idx == 2:
		genome = parts[0]

	# Get chromosome
	chr_name = parts[start_end_idx - 1]

	# Get summit offset
	summit = None
	if start_end_idx + 1 < len(parts):
		summit = parts[start_end_idx + 1]
		# Check that summit offset is an int
		try:
			summit = int(summit)
		except Exception as e:
			print(f"Warning (_get_range): summit offset is not an integer. Got '{name}'")
			return None, None, None, None, None

	return genome, chr_name, start, end, summit



def _check_parallel_files(in_files, out_files):
	"""Check that inputs and outputs have the right numbers of unique files.
	This is intended for functions that take N input files and produce N output files.
	In particular:
	- Check that number of input files equals number of output files
	- Check that output filenames are unique
	"""
	if len(in_files) != len(out_files):
		raise ValueError(f"Number of input files ({len(in_files)}) differs from number of output files ({len(out_files)})")

	if len(set(out_files)) != len(out_files):
		raise ValueError(f"Output filenames are not unique")

def _make_parent_dir(path):
	"""Make parent directory of path"""
	os.makedirs(os.path.abspath(os.path.dirname(path)), exist_ok=True)

def _drop_bed_duplicates(df):
	"""Drop duplicate peaks, assuming a bed/narrowPeak format DataFrame.
	Duplicate key is (chr, start, end), and summit offset if it exists.

	Args:
		df (pd.DataFrame): Peaks, assumed to be in bed/narrowPeak format.
	Returns:
		df (pd.DataFrame): De-duplicated
	"""
	# Duplicate key is (chr, start, end), and summit offset if it exists
	duplicate_key = [0, 1, 2]
	if 9 in df.columns:
		duplicate_key += [9]
	# If there are duplicates, keep only the first one
	return df.sort_values(by=list(df.columns)).drop_duplicates(subset=duplicate_key, keep='first')

def filter_length(in_files, out_files, conditions, combination='and'):
	"""Keep peaks only if their length matches the given conditions.
	Example usage:

	filter_length(..., ..., ["<= 1000"]):
		Keep peaks whose length is <= 1000
	filter_length(..., ..., ["> 2000"]):
		Keep peaks whose length is > 2000
	filter_length(..., ..., [">= 500", "<= 1000"])
		Keep peaks whose length is >= 500 and <= 1000
	filter_length(..., ..., ["< 500", "> 1000"], combination='or')
		Keep peaks whose length is < 500 or > 1000
	"""
	# Check args
	_check_parallel_files(in_files, out_files)
	if combination not in ['and', 'or']:
		raise ValueError(f"Invalid combination '{combination}', valid options are ['and', 'or']")
	if not conditions:
		raise ValueError(f"Missing required parameter --conditions")
	if not isinstance(conditions, list):
		raise ValueError(f"Invalid type for conditions, {type(conditions)}, must be list of str")

	# Get (comparator, value) for each condition
	conditions = [_parse_condition_str(cond) for cond in conditions]
	# Convert comparator strings to functions
	operator_mapping = {"<": operator.lt, ">": operator.gt, "<=": operator.le, ">=": operator.ge}
	conditions = [(operator_mapping[c], v) for (c, v) in conditions]

	messages = []
	for in_file, out_file in zip(in_files, out_files):
		df = pd.read_csv(in_file, delim_whitespace=True, header=None, na_values=['.'])
		orig_len = len(df)

		# Calculate length as start - end
		length = df[2] - df[1]
		# Get rows that match each condition
		match_results = []
		for op, value in conditions:
			match = op(length, value)
			match_results.append(match)
		# Get rows that match all (or any) conditions
		match = pd.concat(match_results, axis=1)
		if combination == 'and':
			match = match.all(axis='columns')
		else:
			match = match.any(axis='columns')
		# Keep only rows that match conditions
		df = df[match]

		messages.append(f"Info (filter_length): Dropped {orig_len - len(df)} peaks, kept {len(df)} of {orig_len}")

		# Make directory and write result .bed or .narrowPeak
		_make_parent_dir(out_file)
		df.to_csv(out_file, index=False, sep="\t", header=None, na_rep='.')

	# Report
	for message in messages:
		print(message)

def _parse_condition_str(condition_str):
	"""E.g. "<= 1000" --> ("<=", 1000)"""

	# Order of comparators is important, to avoid parsing "<=" as "<" "="
	comparators = ["<=", ">=", "<", ">"]
	bad_format_msg = f"Condition string '{condition_str}' has invalid format. Use e.g. '<= 1000'."

	# Remove initial and terminal whitespace
	condition_str = condition_str.strip()

	# Get comparator and value
	valid = False
	e = None
	# Try all possible comparators until we find one that works
	for c in comparators:
		if condition_str.startswith(c):
			comparator = c
			parts = condition_str.split(c)
			if len(parts) == 2:
				try:
					value = int(parts[1])
					valid = True
					break
				except Exception as e:
					bad_format_msg = bad_format_msg + f" {e}"
					raise ValueError(bad_format_msg)

	# Handle errors
	if not valid:
		raise ValueError(bad_format_msg)

	return c, value

def _test_parse_condition_str():
	test_cases = [
		("<= 1000", ("<=", 1000)),
		(">= 1000", (">=", 1000)),
		("< 1000", ("<", 1000)),
		("> 1000", (">", 1000)),
		(" <= 1000  ", ("<=", 1000)),
		(">1000", (">", 1000)),
		("dlfjdk", None),
		("<=1000dlf", None),
		("<=<=1000", None),
		("1000<=", None),
		("<= 1000 1000", None),
		("< = 1000", None)
	]
	for condition_str, expected in test_cases:
		try:
			c, value = _parse_condition_str(condition_str)
			assert (c, value) == expected
		except Exception as e:
			expected is None

def filter_chromosome(in_files, out_files, drop_non_chromosomal=False, drop_values=None):
	"""Filter peaks by chromosome name. Can be used to:
	- Drop peaks that are not on chromosomes (drop_non_chromosomal)
	- Drop peaks whose chromosome name is in a given list (drop_values)

	Args:
		in_files (list of str): can be .bed, .narrowPeak, .bed.gz, or .narrowPeak.gz
		out_files (list of str): can be .bed, .narrowPeak, .bed.gz, or .narrowPeak.gz
		drop_non_chromosomal (bool): if True, then drop rows where chr value is not chrX/chrY/chr<numeric>
		drop_values (list of str): additional chr values to drop
	"""
	# Check args
	_check_parallel_files(in_files, out_files)
	if drop_values is not None:
		if not isinstance(drop_values, list):
			raise ValueError(f"drop_values must be a list, got {type(drop_values)}, {drop_values}")

	messages = []
	for in_file, out_file in zip(in_files, out_files):
		df = pd.read_csv(in_file, delim_whitespace=True, header=None, na_values=['.'])
		orig_len = len(df)

		if drop_non_chromosomal:
			valid = [_is_chromosomal(chr_str) for chr_str in df[0]]
			df = df[valid]

		if drop_values is not None:
			valid = ~df[0].isin(drop_values)
			df = df[valid]

		messages.append(f"Info (filter_chromosome): Dropped {orig_len - len(df)} peaks, kept {len(df)} of {orig_len}")

		# Make directory and write result .bed or .narrowPeak
		_make_parent_dir(out_file)
		df.to_csv(out_file, index=False, sep="\t", header=None, na_rep='.')

	# Report
	for message in messages:
		print(message)

def _is_chromosomal(s):
	"""True if s is 'chrX', 'chrY', or 'chr<numeric>', False otherwise.
	NOTE this should be expanded to include other chromosome naming conventions."""
	if s in ['chrX', 'chrY']:
		return True
	parts = s.split('chr')
	if (len(parts) == 2) and (parts[0] == '') and (parts[1].isnumeric()):
		return True
	return False

def split(in_files, split_ref_species, out_dir,
			split_ref_dir=None, do_get_original=False, create_folds=None):
	"""Split bed files into train, validation, and test sets by Synteny method.

	Args:
		in_files (list of str): Input bed files with regions to split.
			The input bed files should be in the coordinates of one of the five
			reference species.
		split_ref_species (str): Reference species to use for split. Options:
			human, Homo_sapiens, hg38
			mouse, Mus_musculus, mm10
			rhesus, Macaca_mulatta, rheMac10
			rheMac8
			rat, Rattus_norvegicus, rn6
		out_dir (str): Path to output dir. Outputs will have the following structure:
			<out_dir>/
				fold1/
					train.bed
					valid.bed
				[fold2/]
				[...]
				[fold5/]
				test/
					test.bed
		split_ref_dir (str): Dir containing reference .bed files for splitting.
			Default is ../data/bed/ and this is the canonical set of splits.
		do_get_original (bool):
			if False (default), then the output files will be in the same coordinates as the input.
			if True, then the output files will be in the coordinates of the original species,
				taken from the Name column (e.g. chr1:100000-100500)
		create_folds (list of (str or int)):
			Folds to create in a cross-fold validation scheme.
			Default is [1], which will create only fold1 and the test set.
	"""
	# Parse args
	if split_ref_dir is None:
		split_ref_dir = os.path.abspath('../data/bed/')

	if create_folds is None:
		create_folds = [1]
	else:
		create_folds = [int(fold) for fold in create_folds]
	if not all(itm in [1, 2, 3, 4, 5] for itm in create_folds):
		raise ValueError(f"Invalid entry for folds (must be 1, 2, 3, 4, 5) got {create_folds}")

	split_ref_species = SPLIT_REFERENCE_SPECIES_MAPPING[split_ref_species.lower()]

	# Union input files
	import tempfile
	with tempfile.NamedTemporaryFile(suffix='.bed') as union:
		union_bt = subtract(in_files, None, out_files=union.name, window=0)

		# Create test set
		test_ref = _split_reference_bed(split_ref_species, 1, 'test', split_ref_dir)
		test_set = os.path.join(out_dir, 'test', 'test.bed')
		_make_parent_dir(test_set)
		# Take anything that overlaps the test reference
		union_bt.intersect(BedTool(test_ref), wa=True).saveas(test_set)
		# If do_get_original == True, then convert back to the coordinates in the original species
		if do_get_original:
			get_original([test_set], [test_set])

		# For each fold:
		for fold in create_folds:

			# Create validation set
			val_ref = _split_reference_bed(split_ref_species, fold, 'valid', split_ref_dir)
			val_set = os.path.join(out_dir, f'fold{fold}', 'valid.bed')
			_make_parent_dir(val_set)
			# Remove anything that overlaps the test reference,
			# then take anything that overlaps the validation reference
			union_bt.window(BedTool(test_ref), v=True, w=0
				).intersect(BedTool(val_ref), wa=True
				).saveas(val_set)
			if do_get_original:
				get_original([val_set], [val_set])

			# Create training set
			train_ref = _split_reference_bed(split_ref_species, fold, 'train', split_ref_dir)
			train_set = os.path.join(out_dir, f'fold{fold}', 'train.bed')
			_make_parent_dir(train_set)
			# Remove anything that overlaps the test reference,
			# then remove anything that overlaps the validation reference,
			# then take what remains
			union_bt.window(BedTool(test_ref), v=True, w=0
				).window(BedTool(val_ref), v=True, w=0
				).saveas(train_set)
			if do_get_original:
				get_original([train_set], [train_set])

def _split_reference_bed(species, fold, partition, split_ref_dir):
	if species not in SPLIT_REFERENCE_SPECIES_MAPPING.values():
		raise ValueError(f"Invalid species: {species}")
	if fold not in [1, 2, 3, 4, 5]:
		raise ValueError(f"Invalid fold: {fold}")
	if partition not in ['train', 'valid', 'test']:
		raise ValueError(f"Invalid partition: {partition}")
	return os.path.join(split_ref_dir, f"{species}.fold{fold}.{partition}.bed")


def get_args():
	import argparse
	parser = argparse.ArgumentParser()
	# Use subparsers so each function has its own set of arguments
	subparsers = parser.add_subparsers()

	# Arguments for expand_peaks
	parser_expand_peaks = subparsers.add_parser('expand_peaks')
	parser_expand_peaks.set_defaults(function='expand_peaks')
	parser_expand_peaks.add_argument('--in_files', '-i', nargs='+')
	parser_expand_peaks.add_argument('--out_files', '-o', nargs='+')
	parser_expand_peaks.add_argument('--length', '-l', type=int)

	# Arguments for subtract
	parser_subtract = subparsers.add_parser('subtract')
	parser_subtract.set_defaults(function='subtract')
	parser_subtract.add_argument('--include_files', nargs='*')
	parser_subtract.add_argument('--exclude_files', nargs='*')
	parser_subtract.add_argument('--out_files', '-o', nargs='+')
	parser_subtract.add_argument('--window', '-w', type=int)

	# Arguments for get_original
	parser_get_original = subparsers.add_parser('get_original')
	parser_get_original.set_defaults(function='get_original')
	parser_get_original.add_argument('--in_files', '-i', nargs='+')
	parser_get_original.add_argument('--out_files', '-o', nargs='+')

	# Arguments for filter_length
	parser_filter_length = subparsers.add_parser('filter_length')
	parser_filter_length.set_defaults(function='filter_length')
	parser_filter_length.add_argument('--in_files', '-i', nargs='+')
	parser_filter_length.add_argument('--out_files', '-o', nargs='+')
	parser_filter_length.add_argument('--conditions', nargs='+')
	parser_filter_length.add_argument('--combination', default='and', choices=['and', 'or'])

	# Arguments for filter_chromosome
	parser_filter_chromosome = subparsers.add_parser('filter_chromosome')
	parser_filter_chromosome.set_defaults(function='filter_chromosome')
	parser_filter_chromosome.add_argument('--in_files', '-i', nargs='+')
	parser_filter_chromosome.add_argument('--out_files', '-o', nargs='+')
	parser_filter_chromosome.add_argument('--drop_non_chromosomal', action='store_true')
	parser_filter_chromosome.add_argument('--drop_values', nargs='+')

	# Arguments for split
	parser_split = subparsers.add_parser('split')
	parser_split.set_defaults(function='split')
	parser_split.add_argument('--in_files', '-i', nargs='+', required=True)
	parser_split.add_argument('--split_ref_species', '-s', required=True)
	parser_split.add_argument('--split_ref_dir', '-r')
	parser_split.add_argument('--get_original', '-g', action='store_true')
	parser_split.add_argument('--out_dir', '-o', required=True)
	parser_split.add_argument('--create_folds', '-f', nargs='+')

	return parser.parse_args()

if __name__ == '__main__':
	main(get_args())
