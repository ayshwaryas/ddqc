# file to make results_dir available in all functions as a global variable


results_dir = ""


def set_results_dir(d):
    global results_dir
    results_dir = d
