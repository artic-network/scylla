# Tests

To run tests:

`nf-test test --profile docker`

## Tags

Tests are tagged with one of the 3 following tags:

`ci` - Tests that will run as part of the CI pipeline

`unit` - Tests for a specific workflows

`e2e` - End to end tests

This is useful to run a subset of tests.

## Test dataset generation

The test FASTQ files for HCID and Spike sequences were generated using the following process:
- Extract the reference sequence (located in `resources/hcid_refs.fa.gz`)
- Randomly subset to 10 reads of 1000bp length from the reference
- Set quality scores artificially to 40 for all bases
- Added 10 additional reads of random nucleotides (ATCG) at 1000bp length

Each generated FASTQ file was validated by:
- Mapping back to its reference using the same minimap2 settings as the Scylla pipeline
- Classification with Kraken2 + PlusPF database

The test dataset should produce:
- 10 HCID reads identified as CCHF (detected by both Kraken2 and mapping)
- 10 reads that are unclassified (the random reads)

To reduce the size, the NCBI taxonomy dump was downsampled to contain only CCHF with the following commands:

```bash
taxonkit --data-dir taxonomy_dir list --ids 3052518 -I "" \
    | taxonkit --data-dir taxonomy_dir filter -E species \
    | taxonkit --data-dir taxonomy_dir lineage -t \
    | cut -f 3 \
    | sed -s 's/;/\n/g' \
    > taxids.txt

echo 1 >> taxids.txt

mkdir subset
grep -w -f <(awk '{print "^"$1}' taxids.txt) taxonomy_dir/nodes.dmp > subset/nodes.dmp
grep -w -f <(awk '{print "^"$1}' taxids.txt) taxonomy_dir/names.dmp > subset/names.dmp

touch subset/delnodes.dmp subset/merged.dmp
```
## Adding tests

When adding new tests please refer to the [nf-test  documentation](https://www.nf-test.com/) as a starting point.

The [nf-core training module](https://training.nextflow.io/2.4.0/side_quests/nf-test/) also provides a comprehensive walkthrough for writing tests.

### Writing a test
1. Identify what type of component you are testing - e.g process, pipeline, workflow
2. Generate a stub using `nf-test generate`
3. Examine and mock the inputs to that workflow
    - How many channels?
    - What do the channels contain? Tuples? A single path?
    - If it contains tuples - what do they contain?
    - Look at process inputs to help with this
4. Assert agaisnt expected outcomes (mock data or snapshot)

In general [snapshot testing](https://www.nf-test.com/docs/assertions/snapshots/) is preferred due to comparison of md5 checksums which ensures no variability in outputs. In cases where tools produce non-deterministic outputs, you can check for strings in files using `.contains()` or worst case that the file exists with `.exists()`. However the `snapshot > contains > exists` hierrarchy should be mantained where possible.

### Adding tests to CI

Tests in this repo are run on every pull request using GitHub Actions. To ensure any new tests run in this CI pipeline, add the `ci` tag to the test e.g `tag "ci"`. Note that tests can have multiple tags.

## Todo list
- [ ] Illumina tests
- [ ] More comprehensive human detection tests