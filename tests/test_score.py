import swga.score


def test_read_set_finder_line():
    line = "1,2,3,4 3500.01234"
    results = swga.score.read_set_finder_line(line)
    assert results == ([1,2,3,4], 3500.01234)


def test_seq_diff():
    seq = [1, 3, 4, 5]
    unsorted_seq = [3, 1, 5, 4]
    diffs = [2, 1, 1]
    assert swga.score.seq_diff(seq) == diffs
    assert swga.score.seq_diff(unsorted_seq) == diffs


def test_get_user_fun():
    spec_str = "swga.score:get_user_fun"
    fun = swga.score.get_user_fun(spec_str)
    assert fun == swga.score.get_user_fun


def test_default_score_set():
    expression = "fg_dist_mean/bg_dist_mean"
    primer_set = [1,2,3,4]
    primer_locs = [1,3,5,7]
    score, namespace = swga.score.default_score_set(
        expression=expression,
        primer_set=primer_set,
        primer_locs=primer_locs,
        max_dist=2,
        bg_dist_mean=3.0)
    assert score == 2.0/3.0
    assert '__builtins__' not in namespace.keys()



        
