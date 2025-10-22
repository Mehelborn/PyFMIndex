#!/usr/bin/env python

import os
import ctypes


class BwtBlockList(ctypes.Union):
    _fields_ = [("as_nucleotide", ctypes.POINTER()), ("as_amino", ctypes.POINTER())]


class SearchRange(ctypes.Structure):
    _fields_ = [("start_ptr", ctypes.c_uint64), ("end_ptr", ctypes.c_uint64)]


class IndexConfiguration(ctypes.Structure):
    _fields_ = [("")]
    # uint8_t suffixArrayCompressionRatio;
    # uint8_t kmerLengthInSeedTable;
    # enum AwFmAlphabetType alphabetType;
    # bool keepSuffixArrayInMemory;
    # bool storeOriginalSequence;


class CompressedSuffixArray(ctypes.Structure):
    _fields_ = [
        ("value_bit_width", ctypes.c_uint8),
        ("values", ctypes.POINTER(ctypes.c_uint8)),
        ("compressed_byte_length", ctypes.c_uint64),
    ]


class Index(ctypes.Structure):
    _fields_ = [
        ("version_number", ctypes.c_uint32),
        ("feature_flags", ctypes.c_uint32),  # for non user-customizable options.
        ("bwt_length", ctypes.c_uint64),
        ("bwt_block_list", BwtBlockList),
        ("prefix_sums", ctypes.POINTER(ctypes.c_uint64)),
        ("search_range", ctypes.POINTER(SearchRange)),
        ("file_handle", ctypes.POINTER(ctypes.c_void_p)),
        ("config", IndexConfiguration),
        ("file_descriptor", ctypes.c_int),
        ("suffix_array_file_offset", ctypes.c_size_t),
        ("sequence_file_offset", ctypes.c_size_t),
        ("fasta_vector", ctypes.POINTER()),  # fasta lib
        ("suffix_array", CompressedSuffixArray),
    ]


so_path = os.path.abspath("lib/AvxWindowFmIndex/build/libawfmindex.so")
lib = ctypes.cdll.LoadLibrary(so_path)

create_index = lib.awFmCreateIndex
create_index.argtypes = [ctypes.POINTER(ctypes.POINTER(Index))]
