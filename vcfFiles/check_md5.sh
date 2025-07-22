#!/bin/bash

manifest="20220804_manifest.txt"

while read -r filename expected_md5; do
  if [[ -f "$filename" ]]; then
    actual_md5=$(md5 -q "$filename")
    if [[ "$actual_md5" == "$expected_md5" ]]; then
      echo "$filename: OK"
    else
      echo "$filename: FAILED (expected $expected_md5, got $actual_md5)"
    fi
  else
    echo "$filename: NOT FOUND"
  fi
done < "$manifest"
