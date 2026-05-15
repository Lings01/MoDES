"""CI smoke tests: verify package imports and basic API surface."""


def test_import_modes():
    import modes
    assert modes is not None
    assert modes.__version__ == "0.1.0"


def test_public_api_imports():
    from modes import EventCandidateBuilder, MoDEData, MoDES, MoDESResult
    assert MoDES is not None
    assert MoDEData is not None
    assert MoDESResult is not None
    assert EventCandidateBuilder is not None
