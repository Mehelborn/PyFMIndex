from dataclasses import dataclass
from enum import IntEnum
import logging
import os

import ctypes
from pyfmindex import _pyfmindex as _pfmi

logger = logging.getLogger(__name__)

class ReturnCode(IntEnum):
    Success = 1
    FileReadOkay = 2
    FileWriteOkay = 3
    GeneralFailure = -1
    UnsupportedVersionError = -2
    AllocationFailure = -3
    NullPtrError = -4
    SuffixArrayCreationFailure = -5
    IllegalPositionError = -6
    NoFileSrcGiven = -7
    NoDatabaseSequenceGiven = -8
    FileFormatError = -9
    FileOpenFail = -10
    FileReadFail = -11
    FileWriteFail = -12
    ErrorDbSequenceNull = -13
    ErrorSuffixArrayNull = -14
    FileAlreadyExists = -15

class IndexConfiguration:
    def __init__(
        self,
        suffix_array_compression_ratio: int,
        kmer_length_in_seed_table: int,
        alphabet_type: int,
        keep_suffix_array_in_memory: bool,
        store_original_sequence: bool,
    ) -> None:
        self._config = _pfmi._IndexConfiguration(
            suffix_array_compression_ratio,
            kmer_length_in_seed_table,
            alphabet_type,
            keep_suffix_array_in_memory,
            store_original_sequence,
        )

    @property
    def suffix_array_compression_ratio(self) -> int:
        return self._config.suffix_array_compression_ratio

    @suffix_array_compression_ratio.setter
    def suffix_array_compression_ratio(self, value: int) -> None:
        self._config.suffix_array_compression_ratio = value

    @property
    def kmer_length_in_seed_table(self) -> int:
        return self._config.kmer_length_in_seed_table

    @kmer_length_in_seed_table.setter
    def kmer_length_in_seed_table(self, value: int) -> None:
        self._config.kmer_length_in_seed_table = value

    @property
    def alphabet_type(self) -> int:
        return self._config.alphabet_type

    @alphabet_type.setter
    def alphabet_type(self, value: int) -> None:
        self._config.alphabet_type = value

    @property
    def keep_suffix_array_in_memory(self) -> bool:
        return self._config.keep_suffix_array_in_memory

    @keep_suffix_array_in_memory.setter
    def keep_suffix_array_in_memory(self, value: bool) -> None:
        self._config.keep_suffix_array_in_memory = value

    @property
    def store_original_sequence(self) -> bool:
        return self._config.store_original_sequence

    @store_original_sequence.setter
    def store_original_sequence(self, value: bool) -> None:
        self._config.store_original_sequence = value


@dataclass
class SearchRange:
    start_ptr: int
    end_ptr: int


@dataclass
class KmerSearchData:
    kmer_string: str
    kmer_length: int
    position_list: list
    count: int
    capacity: int


class Index:
    _index = None

    def __init__(
        self,
        config: IndexConfiguration,
        file_path: str,
        sequence: str = None,
        fasta_path: str = None,
    ) -> None:
        if not all(( config.alphabet_type, config.keep_suffix_array_in_memory,
                config.kmer_length_in_seed_table, config.store_original_sequence,
                config.suffix_array_compression_ratio)):  # fmt: skip
            raise Exception("Index configuration is not fully initialized.")

        if os.path.exists(file_path):
            file_path_bytes = file_path.encode()
        else:
            raise FileNotFoundError(f"{file_path} not found.")

        index_ptr = ctypes.POINTER(_pfmi._Index)()
        if sequence and not fasta_path:
            sequence_bytes = sequence.encode()
            sequence_bytes_length = len(sequence_bytes)
            sequence_array = (ctypes.c_uint8 * sequence_bytes_length).from_buffer_copy(
                sequence_bytes
            )

            return_code: int = _pfmi._create_index(
                ctypes.byref(index_ptr),
                ctypes.byref(config._config),
                sequence_array,
                len(sequence_bytes),
                file_path_bytes,
            )
            if return_code != ReturnCode.FileWriteOkay:
                raise Exception(f"Failed to create index: {ReturnCode(return_code).name}.")

        elif fasta_path:
            if os.path.exists(fasta_path):
                fasta_path_bytes = fasta_path.encode()
            else:
                raise FileNotFoundError(f"{fasta_path}")

            return_code: int = _pfmi._create_index_from_fasta(
                ctypes.byref(index_ptr),
                ctypes.byref(config._config),
                fasta_path_bytes,
                file_path_bytes,
            )
            if return_code != ReturnCode.FileWriteOkay:
                raise Exception(f"Failed to create index from fasta: {ReturnCode(return_code).name}.")
        else:
            raise Exception("No data for index building.")

        self._index = index_ptr

    def find_search_range_for_string(self, kmer: str) -> SearchRange | None:
        kmer_bytes = kmer.encode()
        kmer_length = len(kmer_bytes)
        if kmer_length == 0:
            raise ValueError("Invalid length")
        search_range: _pfmi._SearchRange = _pfmi._find_search_range_for_string(
            self._index, kmer_bytes, kmer_length
        )
        if search_range.start_ptr >= search_range.end_ptr:
            return None
        return SearchRange(search_range.start_ptr, search_range.end_ptr)

    def read_sequence_from_file(self, start: int, segment_length: int):
        if start < 0:
            raise ValueError("The start position must be more or equal to 0")
        if segment_length <= 0:
            return ""
        buffer_size = segment_length + 1
        buffer = ctypes.create_string_buffer(buffer_size)
        return_code: int = _pfmi._read_sequence_from_file(
            self._index, start, segment_length, buffer
        )
        if return_code == ReturnCode.FileReadFail:
            raise IOError(
                f"ERROR: {ReturnCode.FileReadFail.name}. Could not read the index file."
            )
        elif return_code == ReturnCode.IllegalPositionError:
            raise ValueError(
                f"ERROR: {ReturnCode.IllegalPositionError.name}. The start position is not less than the end position."
            )
        elif return_code == ReturnCode.UnsupportedVersionError:
            raise Exception(
                f"ERROR: {ReturnCode.UnsupportedVersionError.name}. The index was configured to not store the original sequence."
            )
        return buffer.value.decode()

    @property
    def version_number(self):
        return self._index.contents.version_number

    @property
    def feature_flags(self):
        return self._index.contents.feature_flags

    @property
    def bwt_length(self):
        return self._index.contents.bwt_length

    @property
    def bwt_block_list(self):
        return self._index.contents.bwt_block_list

    @property
    def prefix_sums(self):
        return self._index.contents.prefix_sums

    @property
    def kmer_seed_table(self):
        return self._index.contents.kmer_seed_table

    @property
    def file_handle(self):
        return self._index.contents.file_handle

    @property
    def config(self):
        return self._index.contents.config

    @property
    def file_descriptor(self):
        return self._index.contents.file_descriptor

    @property
    def suffix_array_file_offset(self):
        return self._index.contents.suffix_array_file_offset

    @property
    def sequence_file_offset(self):
        return self._index.contents.sequence_file_offset

    @property
    def fasta_vector(self):
        return self._index.contents.fasta_vector

    @property
    def suffix_array(self):
        return self._index.contents.suffix_array

    def __del__(self):
        if self._index is not None:
            _pfmi._dealloc_index(self._index)


def read_index_from_file(file_path: str, keep_suffix_array_in_memory: bool = False):
    if os.path.exists(file_path):
        file_path_bytes = file_path.encode()
    else:
        raise FileNotFoundError(f"{file_path} not found.")
    index_ptr = POINTER(_pfmi._Index)()

    return_code: int = _pfmi._read_index_from_file(
        byref(index_ptr),
        file_path_bytes,
        keep_suffix_array_in_memory,
    )
    if return_code != ReturnCode.FileReadOkay:
        raise Exception(f"Failed to read index from file: {ReturnCode(return_code).name}")
    return index_ptr


class KmerSearchList:
    def __init__(self, capacity: int) -> None:
        if capacity <= 0:
            raise ValueError("Invalid capacity")
        self._kmer_search_list = _pfmi._create_kmer_search_list(capacity)

    def fill_out_list(self, kmers: list[str], kmers_length: list[int]):
        num_kmers = len(kmers)
        kmer_search_data = self._kmer_search_list.contents.kmer_search_data
        for i in range(num_kmers):
            kmer_search_data[i].kmer_string = kmers[i].encode()
            kmer_search_data[i].kmer_length = kmers_length[i]
        self._kmer_search_list.contents.count = num_kmers

    @property
    def capacity(self):
        return self._kmer_search_list.contents.capacity

    @property
    def count(self):
        return self._kmer_search_list.contents.count

    @property
    def kmer_search_data(self):
        return self._kmer_search_list.contents.kmer_search_data

    def __del__(self):
        if self._kmer_search_list:
            _pfmi._dealloc_kmer_search_list(self._kmer_search_list)


def parallel_search_locate(
    index: Index, kmer_search_list: KmerSearchList, num_threads: int = 4
):
    return_code = _pfmi._parallel_search_locate(
        index._index, kmer_search_list._kmer_search_list, num_threads
    )
    if return_code != ReturnCode.Success:
        raise Exception(f"Failed to find the k-mers locations: {ReturnCode(return_code).name}")

def parallel_search_count(
    index: Index, kmer_search_list: KmerSearchList, num_threads: int = 4
):
    _pfmi._parallel_search_count(
        index._index, kmer_search_list._kmer_search_list, num_threads
    ) # awFmParallelSearchCount returns void and fills KmerSearchData
