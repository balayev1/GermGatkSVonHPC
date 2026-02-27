#!/usr/bin/env python

import argparse
import json
import re
from pathlib import Path


def parse_pairs(items):
    pairs = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"Invalid key=value pair: {item}")
        key, value = item.split("=", 1)
        pairs[key] = value
    return pairs


def find_placeholders(value, path=""):
    errors = []
    if isinstance(value, dict):
        for key, item in value.items():
            child_path = f"{path}.{key}" if path else str(key)
            errors.extend(find_placeholders(item, child_path))
    elif isinstance(value, list):
        for idx, item in enumerate(value):
            child_path = f"{path}[{idx}]"
            errors.extend(find_placeholders(item, child_path))
    elif isinstance(value, str) and re.fullmatch(r"__[^_].*__", value):
        errors.append((path, value))
    return errors


def main():
    parser = argparse.ArgumentParser(description="Render JSON template with static/dynamic overrides.")
    parser.add_argument("--template", required=True, help="Path to JSON template")
    parser.add_argument("--out", required=True, help="Output JSON file path")
    parser.add_argument("--static-json", default="{}", help="JSON object merged into template")
    parser.add_argument("--merge-json-file", action="append", default=[], help="Path to JSON object to merge into template")
    parser.add_argument("--set", action="append", default=[], help="Set string key=value")
    parser.add_argument("--set-json", action="append", default=[], help="Set JSON key=<json-value>")
    args = parser.parse_args()

    template_data = json.loads(Path(args.template).read_text())
    static_data = json.loads(args.static_json)
    if not isinstance(static_data, dict):
        raise ValueError("--static-json must be a JSON object")

    rendered = dict(template_data)
    rendered.update(static_data)

    for merge_file in args.merge_json_file:
        merge_data = json.loads(Path(merge_file).read_text())
        if not isinstance(merge_data, dict):
            raise ValueError(f"--merge-json-file must contain a JSON object: {merge_file}")
        rendered.update(merge_data)

    rendered.update(parse_pairs(args.set))

    for key, raw_value in parse_pairs(args.set_json).items():
        rendered[key] = json.loads(raw_value)

    unresolved = find_placeholders(rendered)
    if unresolved:
        details = ", ".join([f"{path}={value}" for path, value in unresolved])
        raise ValueError(f"Unresolved template placeholders found: {details}")

    Path(args.out).write_text(json.dumps(rendered, indent=2) + "\n")


if __name__ == "__main__":
    main()
