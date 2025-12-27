def strip_string(string):
    nlen = len(string)

    k2 = 0
    for idx in range(nlen, 0, -1):
        if string[idx - 1] != " ":
            k2 = idx
            break

    k1 = 1
    for idx in range(1, k2 + 1):
        if string[idx - 1] != " ":
            k1 = idx
            break

    ns = k2 - k1 + 1
    if ns == 0:
        return string, 0

    stripped = string[k1 - 1 : k2]
    stripped = stripped + " " * (nlen - len(stripped))
    return stripped, ns
