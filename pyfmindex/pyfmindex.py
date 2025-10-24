import logging

from ctypes import POINTER, byref, c_uint8
from pyfmindex import _pyfmindex as _pfd

logger = logging.getLogger(__name__)


class IndexConfiguration:
    def __init__(
        self,
        suffix_array_compression_ratio: int,
        kmer_length_in_seed_table: int,
        alphabet_type: int,
        keep_suffix_array_in_memory: bool,
        store_original_sequence: bool,
    ) -> None:
        self._config = _pfd._IndexConfiguration(
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


class Index:
    _index = None

    def __init__(
        self,
        config: IndexConfiguration,
        file_path: str,
        sequence: str = None,
        fasta_path: str = None,
    ) -> None:
        # TODO: check if config fully initialized

        # TODO: check if file_path exists
        file_path_bytes = file_path.encode()

        index_ptr = POINTER(_pfd._Index)()
        if not fasta_path:
            sequence_bytes = sequence.encode()
            sequence_bytes_length = len(sequence_bytes)
            sequence_array = (c_uint8 * sequence_bytes_length).from_buffer_copy(
                sequence_bytes
            )

            return_code: int = _pfd._create_index(
                byref(index_ptr),
                byref(config._config),
                sequence_array,
                len(sequence_bytes),
                file_path_bytes,
            )
            # TODO: check return code

        else:
            # TODO: check if fasta path exists
            return_code: int = _pfd._create_index_from_fasta(
                byref(index_ptr),
                byref(config._config),
                fasta_path.encode(),
                file_path_bytes,
            )
            # TODO: check return code
        self._index = index_ptr

    def write_to_file(self, file_path: str):
        # TODO: check if have a sequence
        # TODO: check if file_path exists
        result = _pfd._write_index_to_file(
            byref(self._index),
            self.sequence.encode(),
            self.sequence_length,
            file_path.encode(),
        )
        # TODO: check return code

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
            _pfd._dealloc_index(self._index)


def read_index_from_file(file_path: str, keep_suffix_array_in_memory: bool = False):
    # TODO: check if path exists
    index_ptr = POINTER(_pfd._Index)()
    return_code: int = _pfd._read_index_from_file(
        byref(index_ptr),
        file_path.encode(),
        keep_suffix_array_in_memory,
    )
    # TODO: check return code
    return index_ptr
