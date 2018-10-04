import subprocess as sp
import re
import os


def find_version():
    version_prefix = "0.0.3"
    try:
        # Try to use git to find current commit.
        path = os.path.dirname(os.path.realpath(__file__))
        p = sp.Popen(
            ["git", "describe", "--always"], cwd=path, stdout=sp.PIPE, stderr=sp.PIPE
        )
        output, error = p.communicate()
        git_describe = output.decode("utf-8").strip()
        git_hash = re.sub("-g[0-9a-f]*$", "", git_describe)

        p = sp.Popen(
            ["git", "log", "-1", "--format=%cd", "--date=iso"],
            cwd=path,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
        )
        output, error = p.communicate()
        git_date = output.decode("utf-8").strip()
        date_time = "_".join(git_date.split(" "))

        __version__ = date_time + "-" + git_hash + "-" + version_prefix
    except BaseException:
        __version__ = version_prefix
    return __version__
