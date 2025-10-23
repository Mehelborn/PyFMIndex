from pyfmindex import IndexConfiguration

SEQUENCE = "ACGTACGTTTACAGT"
SUFFIX_ARRAY_COMPRESSION_RATIO=8
KMER_LENGTHIN_SEED_TABLE=12
ALPHABET_TYPE=2


def test_index_config_creation():
    config = IndexConfiguration(SUFFIX_ARRAY_COMPRESSION_RATIO, KMER_LENGTHIN_SEED_TABLE, ALPHABET_TYPE, True, True)
    assert config.suffix_array_compression_ratio
    assert config.kmer_length_in_seed_table
    assert config.alphabet_type
    assert config.keep_suffix_array_in_memory
    assert config.store_original_sequence
