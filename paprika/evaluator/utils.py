import pkg_resources


def get_benchmarks():
    """
    Determine the installed ``taproom`` benchmarks.
    """
    installed_benchmarks = {}

    for entry_point in pkg_resources.iter_entry_points(group="taproom.benchmarks"):
        installed_benchmarks[entry_point.name] = entry_point.load()

    return installed_benchmarks
