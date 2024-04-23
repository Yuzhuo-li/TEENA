#!/bin/bash

# 检查输入参数
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <bed_file>"
    exit 1
fi

# 输入的BED文件
bed_file=$1

# 检查文件是否存在
if [ ! -f "$bed_file" ]; then
    echo "Error: File '$bed_file' not found!"
    exit 1
fi

# 创建一个新的输出文件名
output_file="${bed_file%.bed}.modified.bed"

# 处理BED文件
awk 'BEGIN {OFS="\t"} {
    if (split($4, arr, ":") == 1 && split($4, arr, "/") == 2) {
        $4 = $4 ":NA"  # 如果第四列有两个部分通过 / 分隔，但没有 :，添加 :NA
    }
    print $0
}' "$bed_file" > "$output_file"

echo "Modified BED file created: $output_file"

