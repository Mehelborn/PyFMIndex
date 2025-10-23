from ctypes import POINTER, byref
from pyfmindex import _pyfmindex as _pfd


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

    @property
    def kmer_length_in_seed_table(self) -> int:
        return self._config.kmer_length_in_seed_table

    @property
    def alphabet_type(self) -> int:
        return self._config.alphabet_type

    @property
    def keep_suffix_array_in_memory(self) -> bool:
        return self._config.keep_suffix_array_in_memory

    @property
    def store_original_sequence(self) -> bool:
        return self._config.store_original_sequence


class Index:
    _index = None

    def __init__(
        self,
        config: IndexConfiguration,
        file_path: str,
        sequence: str = None,
        fasta_path: str = None,
    ) -> None:
        self.config = config  # TODO: check if fully initialized
        self.file_path = file_path  # TODO: check if file_path exists

        index = _pfd._Index()
        if not fasta_path:
            self.sequence = sequence
            self.sequence_length = len(sequence)  # NOTE: check if empty?
            return_code: int = _pfd._create_index(
                byref(POINTER(index)()),
                byref(config),
                sequence.encode(),
                self.sequence_length,
                file_path.encode(),
            )
            # TODO: check return code

        else:
            # TODO: check if fasta path exists
            return_code: int = _pfd._create_index_from_fasta(
                byref(POINTER(index)()),
                byref(config),
                fasta_path.encode(),
                file_path.encode(),
            )
            # TODO: check return code
        self._index = index

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
        return self._index.version_number

    @property
    def feature_flags(self):
        return self._index.feature_flags

    @property
    def bwt_length(self):
        return self._index.bwt_length

    @property
    def bwt_block_list(self):
        return self._index.bwt_block_list

    @property
    def prefix_sums(self):
        return self._index.prefix_sums

    @property
    def kmer_seed_table(self):
        return self._index.kmer_seed_table

    @property
    def file_handle(self):
        return self._index.file_handle

    @property
    def config(self):
        return self._index.config

    @property
    def file_descriptor(self):
        return self._index.file_descriptor

    @property
    def suffix_array_file_offset(self):
        return self._index.suffix_array_file_offset

    @property
    def sequence_file_offset(self):
        return self._index.sequence_file_offset

    @property
    def fasta_vector(self):
        return self._index.fasta_vector

    @property
    def suffix_array(self):
        return self._index.suffix_array

    def __del__(self):
        if self._index is not None:
            _pfd._dealloc_index(byref(self._index))


def read_index_from_file(
    file_path: str, keep_suffix_array_in_memory: bool = True
) -> Index:
    # TODO: check if path exists
    index = _pfd._Index()  # TODO: convert to python Index
    return_code: int = _pfd._read_index_from_file(
        byref(POINTER(index)()),
        file_path.encode(),
        keep_suffix_array_in_memory,
    )
    # TODO: check return code
    return index
