import re
import math


def disect_line(line):
    raw = line.strip("\n")
    ri = raw[:6].strip()
    raw_n = raw[6:]
    matched = re.search(
        r"(\d+)\s*([A-Za-z0-9]+\'?)\s+([A-Za-z0-9]{1,4})\s*([A-Z0-9])\s*(\w*\d+\w*)\s+(-?\d+.?\d+)\s*(-?\d+.?\d+)\s*(-?\d+.?\d+)\s+(\d+.?\d+)\s*(\d+.?\d+)\s+(\w+)",
        raw_n,
    )
    # print(matched.group(1))
    if matched and ri in ["ATOM", "HETATM"]:
        # print(raw)
        p = [
            ri,
            matched.group(1),
            matched.group(2),
            matched.group(3),
            matched.group(4),
            matched.group(5),
            float(matched.group(6)),
            float(matched.group(7)),
            float(matched.group(8)),
            float(matched.group(9)),
            matched.group(10),
            matched.group(11),
        ]
        # print(p)
        return p


def read_pdb(file):
    if file.endswith(".pdb"):
        with open(file) as fh:
            lines = []
            parsed = []
            for line in fh.readlines():
                raw = line.strip("\n")
                ri = raw[:6].strip()
                raw_n = raw[6:]
                matched = re.search(
                    r"(\d+)\s*([A-Za-z0-9]+\'?)\s+([A-Za-z0-9]{1,4})\s*([A-Z0-9])\s*(\w*\d+\w*)\s+(-?\d+.?\d+)\s*(-?\d+.?\d+)\s*(-?\d+.?\d+)\s+(\d+.?\d+)\s*(\d+.?\d+)\s+(\w+)",
                    raw_n,
                )
                # print(matched.group(1))
                if matched and ri in ["ATOM", "HETATM"]:
                    # print(raw)
                    p = [
                        ri,
                        matched.group(1),
                        matched.group(2),
                        matched.group(3),
                        matched.group(4),
                        matched.group(5),
                        float(matched.group(6)),
                        float(matched.group(7)),
                        float(matched.group(8)),
                        float(matched.group(9)),
                        matched.group(10),
                        matched.group(11),
                    ]
                    # print(p)
                    parsed.append(p)
                    lines.append(raw)
        # print(len(parsed), type(parsed))
        return parsed, lines


def distance(p1, p2):
    x = p1[6] - p2[6]
    y = p1[7] - p2[7]
    z = p1[8] - p2[8]
    d = math.sqrt(x * x + y * y + z * z)
    return d


# def main():
#     for file in os.listdir('/Users/kanavkalra/Desktop/all_pdb/'):
#      #for file in ['3u5h.pdb']:
#         if file.endswith('.pdb'):
#             print(file)
#             parsed, lines = read_pdb('/home/kalrak/all_pdb/' + file)
#             print(parsed)


# if __name__ == '__main__':
#     main()
