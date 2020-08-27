#!/usr/bin/env bash

# Download PDB collection
echo "$(date) :: Downloading PDB files"
rsync -rlpt -z --delete --port=33444 --info=progress2 \
rsync.rcsb.org::ftp_data/structures/divided/pdb/ ./pdb

echo "$( date ) :: Downloading mmCIF files"
rsync -rlpt -z --delete --port=33444 --info=progress2 \
rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ ./mmCIF

echo "$( date ) :: Downloading XML files"  # for validation
rsync -rlpt -z --delete --port=33444 --info=progress2 \
rsync.rcsb.org::ftp_data/structures/divided/XML/ ./XML

