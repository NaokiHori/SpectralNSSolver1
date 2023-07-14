import os
import sys
import glob
import enum

def get_filenames(root):
    results = glob.glob(f"{root}/**", recursive=True)
    retvals = list()
    for result in results:
        if result.endswith(".c") or result.endswith(".h"):
            retvals.append(result)
    return retvals

class NdimsType(enum.Enum):
    OUT   = enum.auto()
    IN_2D = enum.auto()
    IN_3D = enum.auto()

def extract(ndims, lines):
    state = NdimsType.OUT
    if_level = 0
    if_level_ndims = 0
    newlines = list()
    for line in lines:
        is_on_ndims_macro = False
        if "#if" in line:
            # found "if", increase nest counter
            if_level += 1
            if " NDIMS" in line:
                if "NDIMS==2" in line.replace(" ", ""):
                    is_on_ndims_macro = True
                    # now in 2D condition
                    state = NdimsType.IN_2D
                    if_level_ndims = if_level
                if "NDIMS==3" in line.replace(" ", ""):
                    is_on_ndims_macro = True
                    # now in 3D condition
                    state = NdimsType.IN_3D
                    if_level_ndims = if_level
        elif "#else" in line:
            # check this "else" is for ndims
            if if_level == if_level_ndims:
                is_on_ndims_macro = True
                # if it is, swap state (3d if now 2d, vice versa)
                if state == NdimsType.IN_2D:
                    state = NdimsType.IN_3D
                elif state == NdimsType.IN_3D:
                    state = NdimsType.IN_2D
                else:
                    print("although we should be inside ndims condition, state is OUT")
                    sys.exit()
        elif "#endif" in line:
            if if_level == if_level_ndims:
                is_on_ndims_macro = True
                state = NdimsType.OUT
            # found "endif", reduce nest counter
            if_level -= 1
        if not is_on_ndims_macro:
            # we do not include macro about ndims
            if ndims == 2 and state != NdimsType.IN_3D:
                newlines.append(line)
            if ndims == 3 and state != NdimsType.IN_2D:
                newlines.append(line)
    return newlines

def modify_comment(lines):
    newlines = list()
    for line in lines:
        if "//" in line and "|" in line:
            line = line.split("|")[0] + "\n"
        newlines.append(line)
    return newlines

def main():
    argv = sys.argv
    # sanitise input
    assert(len(argv) == 2)
    ndims = int(argv[1])
    # target scripts
    fnames = \
          get_filenames("src") \
        + get_filenames("include") \
        + get_filenames("initial_conditions/src") \
        + get_filenames("initial_conditions/include")
    for fname in fnames:
        with open(fname, "r") as f:
            lines = f.readlines()
        lines = extract(ndims, lines)
        lines = modify_comment(lines)
        if 0 == len(lines):
            os.system(f"rm {fname}")
            continue
        if "extern int dummy" in "".join(lines):
            os.system(f"rm {fname}")
            continue
        with open(fname, "w") as f:
            for line in lines:
                f.write(line)

if __name__ == "__main__":
    main()

