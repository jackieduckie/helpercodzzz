#!/bin/bash

vcf_path=\"$1\"
normalName=\"$2\"
genomeVersion=\"hg19\"
module load rc-tools

expr=$(cat <<EOF
  with import <bionix> {} ; (callBionixE (/. + <nixpkgs> + "/../gridss-somatic-filter.nix")
    { normalName = $normalName; genomeVersion = $genomeVersion; }
    $vcf_path).overrideAttrs (_:
    { PPN = 1; MEMORY = "10G"; WALLTIME = "5:00:00"; })
EOF
)

exec nix build -I /stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/stafford/src/nix \
    --impure --expr "$expr" &&

# after build succeeds, run the cmd below and specify the result path
result_path=$3
nix-chroot cp -L ./result $result_path && rm ./result && nix store gc
