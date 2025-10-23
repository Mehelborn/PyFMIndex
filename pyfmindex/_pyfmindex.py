# ruff: noqa F403, F405

import os
from ctypes import *

so_path = os.path.abspath("lib/AvxWindowFmIndex/build/libawfmindex.so")
lib = cdll.LoadLibrary(so_path)


class ReturnCode(c_int):
    Success = 1
    FileReadOkay = 2
    FileWriteOkay = 3
    GeneralFailure = -1
    FileError = -2
    AllocationFailure = -3
    SuffixArrayFailure = -4
    OutOfMemory = -5


class AlphabetType(c_int):
    Amino = 1
    DNA = 2
    RNA = 3


class SimdVec256(Structure):
    _fields_ = [
        ("low_vec", c_uint16),
        ("high_vec", c_uint16),
    ]


class AminoBlock(Structure):
    _fields_ = [
        ("letter_bit_vectors", POINTER(SimdVec256)),
        ("base_occurrences", c_uint64),
    ]


class NucleotideBlock(Structure):
    _fields_ = [
        ("letter_bit_vectors", POINTER(SimdVec256)),
        ("base_occurrences", c_uint64),
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
        ("alphabet_type", c_uint8),
        ("keep_suffix_array_in_memory", c_bool),
        ("store_original_sequence", c_bool),
    ]


class CompressedSuffixArray(Structure):
    _fields_ = [
        ("value_bit_width", c_uint8),
        ("values", POINTER(c_uint8)),
        ("compressed_byte_length", c_uint64),
    ]


class SearchRange(Structure):
    _fields_ = [("start_ptr", c_uint64), ("end_ptr", c_uint64)]


class FastaVectorMetadata(Structure):
    _fields_ = [
        ("headerEndPosition", c_uint32),
        ("sequenceEndPosition", c_uint32),
    ]


class FastaVectorMetadataVector(Structure):
    _fields_ = [
        ("data", POINTER(FastaVectorMetadata)),
        ("capacity", c_uint32),
        ("count", c_uint32),
    ]


class FastaVectorString(Structure):
    _fields_ = [
        ("char_data", c_char_p),
        ("capacity", c_uint32),
        ("count", c_uint32),
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
        ("kmer_seed_table", POINTER(SearchRange)),
        ("file_handle", POINTER(c_void_p)),
        ("config", _IndexConfiguration),
        ("file_descriptor", c_int),
        ("suffix_array_file_offset", c_size_t),
        ("sequence_file_offset", c_size_t),
        ("fasta_vector", POINTER(FastaVector)),
        ("suffix_array", CompressedSuffixArray),
    ]


class KmerSearchData(Structure):
    _fields_ = [
        ("kmer_string", POINTER(c_char)),
        ("kmer_length", c_uint64),
        ("position_list", POINTER(c_uint64)),
        ("count", c_uint32),
        ("capacity", c_uint32),
    ]


class KmerSearchList(Structure):
    _fields_ = [
        ("capacity", c_size_t),
        ("count", c_size_t),
        ("kmer_search_data", POINTER(KmerSearchData)),
    ]


def _is_pointer_to_pointer(obj):
    if isinstance(obj, _Pointer):
        return isinstance(obj._type_, type) and issubclass(obj._type_, _Pointer)
    return False


_create_index = lib.awFmCreateIndex
_create_index.argtypes = [
    POINTER(POINTER(_Index)),
    POINTER(_IndexConfiguration),
    POINTER(c_uint8),
    c_size_t,
    POINTER(c_char),
]


_create_index_from_fasta = lib.awFmCreateIndexFromFasta
_create_index_from_fasta.argtypes = [
    POINTER(POINTER(_Index)),
    POINTER(_IndexConfiguration),
    POINTER(c_char),
    POINTER(c_char),
]

_dealloc_index = lib.awFmDeallocIndex
_dealloc_index.argtypes = [POINTER(_Index)]
_dealloc_index.restype = None

_write_index_to_file = lib.awFmWriteIndexToFile
_write_index_to_file.argtypes = [
    POINTER(_Index),
    POINTER(c_uint8),
    c_uint64,
    POINTER(c_char),
]

_read_index_from_file = lib.awFmReadIndexFromFile
_read_index_from_file.argtypes = [
    POINTER(POINTER(_Index)),
    c_char_p,
    c_bool,
]


find_search_range_for_string = lib.awFmFindSearchRangeForString
find_search_range_for_string.argtypes = [
    POINTER(POINTER(_Index)),
    POINTER(c_char),
    c_size_t,
]
find_search_range_for_string.restype = SearchRange

create_kmer_search_list = lib.awFmCreateKmerSearchList
create_kmer_search_list.argtypes = [c_size_t]
create_kmer_search_list.restype = KmerSearchList

dealloc_kmer_search_list = lib.awFmDeallocKmerSearchList
dealloc_kmer_search_list.argtypes = [POINTER(KmerSearchList)]
dealloc_kmer_search_list.restype = None

parallel_search_locate = lib.awFmParallelSearchLocate
parallel_search_locate.argtypes = [
    POINTER(_Index),
    POINTER(KmerSearchList),
    c_uint32,
]

parallel_search_count = lib.awFmParallelSearchCount
parallel_search_count.argtypes = [
    POINTER(_Index),
    POINTER(KmerSearchList),
    c_uint32,
]


read_sequence_from_file = lib.awFmReadSequenceFromFile
read_sequence_from_file.argtypes = [
    POINTER(_Index),
    c_size_t,
    c_size_t,
    POINTER(c_char),
]

create_initial_query_range = lib.awFmCreateInitialQueryRange
create_initial_query_range.argtypes = [
    POINTER(_Index),
    POINTER(c_char),
    c_uint64,
]
create_initial_query_range.restype = SearchRange

create_initial_query_range_from_char = lib.awFmCreateInitialQueryRangeFromChar
create_initial_query_range_from_char.argtypes = [
    POINTER(_Index),
    c_char,
]
create_initial_query_range_from_char.restype = SearchRange


nucleotide_iterative_step_backward_search = (
    lib.awFmNucleotideIterativeStepBackwardSearch
)
nucleotide_iterative_step_backward_search.argtypes = [
    POINTER(_Index),
    POINTER(SearchRange),
    c_uint8,
]
nucleotide_iterative_step_backward_search.restype = None

amino_iterative_step_backward_search = lib.awFmAminoIterativeStepBackwardSearch
amino_iterative_step_backward_search.argtypes = [
    POINTER(_Index),
    POINTER(SearchRange),
    c_uint8,
]
amino_iterative_step_backward_search.restype = None


find_database_hit_positions = lib.awFmFindDatabaseHitPositions
find_database_hit_positions.argtypes = [
    POINTER(_Index),
    POINTER(SearchRange),
    POINTER(c_int),
]

find_database_hit_position_single = lib.awFmFindDatabaseHitPositionSingle
find_database_hit_position_single.argtypes = [
    POINTER(_Index),
    c_uint64,
    POINTER(c_int),
]
find_database_hit_position_single.restype = POINTER(c_uint64)


get_local_sequence_position_from_index_position = (
    lib.awFmGetLocalSequencePositionFromIndexPosition
)
get_local_sequence_position_from_index_position.argtypes = [
    POINTER(_Index),
    c_size_t,
    POINTER(c_size_t),
    POINTER(c_size_t),
]

nucleotide_backtrace_return_previous_letter_index = (
    lib.awFmNucleotideBacktraceReturnPreviousLetterIndex
)
nucleotide_backtrace_return_previous_letter_index.argtypes = [
    POINTER(_Index),
    POINTER(c_uint64),
]

amino_backtrace_return_previous_letter_index = (
    lib.awFmAminoBacktraceReturnPreviousLetterIndex
)
amino_backtrace_return_previous_letter_index.argtypes = [
    POINTER(_Index),
    POINTER(c_uint64),
]

get_header_string_from_sequence_number = lib.awFmGetHeaderStringFromSequenceNumber
get_header_string_from_sequence_number.argtypes = [
    POINTER(_Index),
    c_size_t,
    POINTER(POINTER(c_char)),
    POINTER(c_size_t),
]

search_range_length = lib.awFmSearchRangeLength
search_range_length.argtypes = [POINTER(SearchRange)]
search_range_length.restype = c_size_t

return_code_is_failure = lib.awFmReturnCodeIsFailure
return_code_is_failure.argtypes = [c_int]
return_code_is_failure.restype = c_bool

return_code_is_success = lib.awFmReturnCodeIsSuccess
return_code_is_success.argtypes = [c_int]
return_code_is_success.restype = c_bool


get_num_sequences = lib.awFmGetNumSequences
get_num_sequences.argtypes = [POINTER(_Index)]
