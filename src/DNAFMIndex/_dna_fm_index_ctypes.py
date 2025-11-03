# ruff: noqa F403, F405

import os
from pathlib import Path
from ctypes import *

AMINO_VECTORS_PER_WINDOW = 5
NUCLEOTIDE_VECTORS_PER_WINDOW = 3
NUCLEOTIDE_CARDINALITY = 4
AMINO_CARDINALITY = 20

ALPHABET_TYPE_AMINO = 1
ALPHABET_TYPE_DNA = 2
ALPHABET_TYPE_RNA = 3

try:
    _lib_dir = Path(__file__).parent / "lib"
except NameError:
    _lib_dir = Path("./src/DNAFMIndex/lib").resolve()

if os.path.exists(_lib_dir):
    _awfmindex = cdll.LoadLibrary(_lib_dir / "libawfmindex.so")
else:
    raise OSError("lib directory not found")


SimdVec256 = c_uint8 * 32


class AminoBlock(Structure):
    _fields_ = [
        ("letter_bit_vectors", SimdVec256 * AMINO_VECTORS_PER_WINDOW),
        (
            "base_occurrences",
            c_uint64 * (AMINO_CARDINALITY + 4),
        ),  # +4 is for sentinel count and 32B padding
    ]


class NucleotideBlock(Structure):
    _fields_ = [
        ("letter_bit_vectors", SimdVec256 * NUCLEOTIDE_VECTORS_PER_WINDOW),
        (
            "base_occurrences",
            c_uint64 * (NUCLEOTIDE_CARDINALITY + 4),
        ),  # +4 is for sentinel count and 32B padding
    ]


class BwtBlockList(Union):
    _fields_ = [
        ("as_nucleotide", POINTER(NucleotideBlock)),
        ("as_amino", POINTER(AminoBlock)),
    ]


class _IndexConfiguration(Structure):
    _fields_ = [
        ("suffix_array_compression_ratio", c_uint8),
        ("kmer_length_in_seed_table", c_uint8),
        ("alphabet_type", c_int),
        ("keep_suffix_array_in_memory", c_bool),
        ("store_original_sequence", c_bool),
    ]


class CompressedSuffixArray(Structure):
    _fields_ = [
        ("value_bit_width", c_uint8),
        ("values", POINTER(c_uint8)),
        ("compressed_byte_length", c_uint64),
    ]


class _SearchRange(Structure):
    _fields_ = [("start_ptr", c_uint64), ("end_ptr", c_uint64)]


class FastaVectorMetadata(Structure):
    _fields_ = [
        ("headerEndPosition", c_size_t),
        ("sequenceEndPosition", c_size_t),
    ]


class FastaVectorMetadataVector(Structure):
    _fields_ = [
        ("data", POINTER(FastaVectorMetadata)),
        ("capacity", c_size_t),
        ("count", c_size_t),
    ]


class FastaVectorString(Structure):
    _fields_ = [
        ("char_data", c_char_p),
        ("capacity", c_size_t),
        ("count", c_size_t),
    ]


class FastaVector(Structure):
    _fields_ = [
        ("sequence", FastaVectorString),
        ("header", FastaVectorString),
        ("metadata", FastaVectorMetadataVector),
    ]


class _Index(Structure):
    _fields_ = [
        ("version_number", c_uint32),
        ("feature_flags", c_uint32),  # for non user-customizable options
        ("bwt_length", c_uint64),
        ("bwt_block_list", BwtBlockList),
        ("prefix_sums", POINTER(c_uint64)),
        ("kmer_seed_table", POINTER(_SearchRange)),
        ("file_handle", POINTER(c_void_p)),
        ("config", _IndexConfiguration),
        ("file_descriptor", c_int),
        ("suffix_array_file_offset", c_size_t),
        ("sequence_file_offset", c_size_t),
        ("fasta_vector", POINTER(FastaVector)),
        ("suffix_array", CompressedSuffixArray),
    ]


class _KmerSearchData(Structure):
    _fields_ = [
        ("kmer_string", c_char_p),
        ("kmer_length", c_uint64),
        ("position_list", POINTER(c_uint64)),
        ("count", c_uint32),
        ("capacity", c_uint32),
    ]


class _KmerSearchList(Structure):
    _fields_ = [
        ("capacity", c_size_t),
        ("count", c_size_t),
        ("kmer_search_data", POINTER(_KmerSearchData)),
    ]


_create_index = _awfmindex.awFmCreateIndex
_create_index.argtypes = [
    POINTER(POINTER(_Index)),
    POINTER(_IndexConfiguration),
    POINTER(c_uint8),
    c_size_t,
    c_char_p,
]


_create_index_from_fasta = _awfmindex.awFmCreateIndexFromFasta
_create_index_from_fasta.argtypes = [
    POINTER(POINTER(_Index)),
    POINTER(_IndexConfiguration),
    c_char_p,
    c_char_p,
]


_dealloc_index = _awfmindex.awFmDeallocIndex
_dealloc_index.argtypes = [POINTER(_Index)]
_dealloc_index.restype = None


_read_index_from_file = _awfmindex.awFmReadIndexFromFile
_read_index_from_file.argtypes = [
    POINTER(POINTER(_Index)),
    c_char_p,
    c_bool,
]


_find_search_range_for_string = _awfmindex.awFmFindSearchRangeForString
_find_search_range_for_string.argtypes = [
    POINTER(_Index),
    c_char_p,
    c_size_t,
]
_find_search_range_for_string.restype = _SearchRange


_create_kmer_search_list = _awfmindex.awFmCreateKmerSearchList
_create_kmer_search_list.argtypes = [c_size_t]
_create_kmer_search_list.restype = POINTER(_KmerSearchList)


_dealloc_kmer_search_list = _awfmindex.awFmDeallocKmerSearchList
_dealloc_kmer_search_list.argtypes = [POINTER(_KmerSearchList)]
_dealloc_kmer_search_list.restype = None


_parallel_search_locate = _awfmindex.awFmParallelSearchLocate
_parallel_search_locate.argtypes = [
    POINTER(_Index),
    POINTER(_KmerSearchList),
    c_uint32,
]


_parallel_search_count = _awfmindex.awFmParallelSearchCount
_parallel_search_count.argtypes = [
    POINTER(_Index),
    POINTER(_KmerSearchList),
    c_uint32,
]
_parallel_search_count.restype = None


_read_sequence_from_file = _awfmindex.awFmReadSequenceFromFile
_read_sequence_from_file.argtypes = [
    POINTER(_Index),
    c_size_t,
    c_size_t,
    c_char_p,
]


create_initial_query_range = _awfmindex.awFmCreateInitialQueryRange
create_initial_query_range.argtypes = [
    POINTER(_Index),
    POINTER(c_char_p),
    c_uint64,
]
create_initial_query_range.restype = _SearchRange


create_initial_query_range_from_char = _awfmindex.awFmCreateInitialQueryRangeFromChar
create_initial_query_range_from_char.argtypes = [
    POINTER(_Index),
    c_char_p,
]
create_initial_query_range_from_char.restype = _SearchRange


nucleotide_iterative_step_backward_search = (
    _awfmindex.awFmNucleotideIterativeStepBackwardSearch
)
nucleotide_iterative_step_backward_search.argtypes = [
    POINTER(_Index),
    POINTER(_SearchRange),
    c_uint8,
]
nucleotide_iterative_step_backward_search.restype = None


amino_iterative_step_backward_search = _awfmindex.awFmAminoIterativeStepBackwardSearch
amino_iterative_step_backward_search.argtypes = [
    POINTER(_Index),
    POINTER(_SearchRange),
    c_uint8,
]
amino_iterative_step_backward_search.restype = None


find_database_hit_positions = _awfmindex.awFmFindDatabaseHitPositions
find_database_hit_positions.argtypes = [
    POINTER(_Index),
    POINTER(_SearchRange),
    POINTER(c_int),
]


find_database_hit_position_single = _awfmindex.awFmFindDatabaseHitPositionSingle
find_database_hit_position_single.argtypes = [
    POINTER(_Index),
    c_uint64,
    POINTER(c_int),
]
find_database_hit_position_single.restype = POINTER(c_uint64)


get_local_sequence_position_from_index_position = (
    _awfmindex.awFmGetLocalSequencePositionFromIndexPosition
)
get_local_sequence_position_from_index_position.argtypes = [
    POINTER(_Index),
    c_size_t,
    POINTER(c_size_t),
    POINTER(c_size_t),
]


nucleotide_backtrace_return_previous_letter_index = (
    _awfmindex.awFmNucleotideBacktraceReturnPreviousLetterIndex
)
nucleotide_backtrace_return_previous_letter_index.argtypes = [
    POINTER(_Index),
    POINTER(c_uint64),
]

amino_backtrace_return_previous_letter_index = (
    _awfmindex.awFmAminoBacktraceReturnPreviousLetterIndex
)
amino_backtrace_return_previous_letter_index.argtypes = [
    POINTER(_Index),
    POINTER(c_uint64),
]


get_header_string_from_sequence_number = (
    _awfmindex.awFmGetHeaderStringFromSequenceNumber
)
get_header_string_from_sequence_number.argtypes = [
    POINTER(_Index),
    c_size_t,
    POINTER(POINTER(c_char_p)),
    POINTER(c_size_t),
]


search_range_length = _awfmindex.awFmSearchRangeLength
search_range_length.argtypes = [POINTER(_SearchRange)]
search_range_length.restype = c_size_t


return_code_is_failure = _awfmindex.awFmReturnCodeIsFailure
return_code_is_failure.argtypes = [c_int]
return_code_is_failure.restype = c_bool


return_code_is_success = _awfmindex.awFmReturnCodeIsSuccess
return_code_is_success.argtypes = [c_int]
return_code_is_success.restype = c_bool


get_num_sequences = _awfmindex.awFmGetNumSequences
get_num_sequences.argtypes = [POINTER(_Index)]
