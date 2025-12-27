import argparse
import json
import sys

from xbl import XFoilState, XBlState, setbl


def _assign_scalar(obj, key, value):
    if hasattr(obj, key):
        setattr(obj, key, value)


def _assign_array(obj, key, value):
    if not hasattr(obj, key):
        return
    arr = getattr(obj, key)
    if isinstance(arr, list):
        if len(arr) != len(value):
            raise ValueError(f"{key}: expected length {len(arr)} got {len(value)}")
        for i, v in enumerate(value):
            arr[i] = v
    else:
        raise ValueError(f"{key}: unsupported array type")


def load_state(path):
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)

    ctx = XFoilState()
    bl = XBlState()

    ctx_data = data.get("ctx", {})
    bl_data = data.get("bl", {})

    for key, value in ctx_data.items():
        if isinstance(value, list):
            _assign_array(ctx, key, value)
        else:
            _assign_scalar(ctx, key, value)

    for key, value in bl_data.items():
        if isinstance(value, list):
            _assign_array(bl, key, value)
        else:
            _assign_scalar(bl, key, value)

    return ctx, bl


def main(argv):
    parser = argparse.ArgumentParser(description="Run Python xbl port on a dumped state.")
    parser.add_argument("state_json", help="Path to JSON dump containing ctx/bl fields.")
    parser.add_argument("--out", help="Optional JSON output path for updated state.")
    args = parser.parse_args(argv)

    ctx, bl = load_state(args.state_json)
    setbl(ctx, bl)

    if args.out:
        with open(args.out, "w", encoding="utf-8") as handle:
            json.dump({"ctx": ctx.__dict__, "bl": bl.__dict__}, handle, indent=2)
    else:
        print(json.dumps({"ctx": ctx.__dict__, "bl": bl.__dict__}, indent=2))


if __name__ == "__main__":
    main(sys.argv[1:])
