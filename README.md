# Evaluation Framework for Data Stream Compression

## Description

***Streaming data*** is one of the primary sources driving today’s big data challenges, due to the immense number of data-generating devices across diverse modern industries. Numerous studies have proposed compression techniques to address the overwhelming volume of such data, with the goal of improving storage efficiency and reducing bandwidth requirements for data transmission. However, such methods are often implemented in different languages and evaluated in disparate environments.

This repository provides a ``unified framework`` for evaluating the performance of compression algorithms in the context of ___univariate___ streaming time series. The framework primarily focuses on comparing lossy functional approximation techniques, due to the fact that they can operate efficiently across arbitrary data domains without requiring explicit training phases.

To ensure a fair evaluation, all algorithms are re-implemented in C++, which is commonly supported by edge devices. Streaming behavior is faithfully simulated by feeding data points to the compressor sequentially. The compressed outputs are serialized in binary format and immediately transmitted to the decompressor for simulating data transfer from sources (e.g., sensors) to sinks (e.g., centralized servers).

## Folder Structure

```text
.
├── bin/              # C++ object and executable files
├── conf/             # Configuration files for each algorithm
├── data/             # Evaluation Datasets
├── include/          # Header files defining [de]compressor and algorithms
├── lib/              # Supporting libraries
├── out/              # Statistical outputs in CSV format
│   ├── compress/     # Compression results
│   └── decompress/   # Decompression results
├── scripts/          # Shell scripts for compilation and execution
└── src/
    ├── C++/          # Algorithm implementations and the main entry point
    └── Python/       # Configuration validation and statistics generation
```

## Environment
- **OS**: `Linux`, preferably an `Ubuntu` based distribution.
- **C++**: The source code is written in `C++11` and compiled using `g++ 11.4.0`.
- **Python**: `numpy` library is required.
- **CPU**: At least ___two cores___ are required due to the use of multithreading for monitoring operations.

## Execution
### Compilation

First, compile all C++ source files by running the following command:

```bash
$ scripts/compile.sh
```

It is worth noting that the project relies heavily on header only libraries. As a result, the compilation process may take a quite amount of time.

### Configuration
Before execution, create a configuration file corresponding to the algorithm wish to run. Configuration files follow the json format below:

```json
{
    "data": Original dataset. Only CSV datasets with two columns (time, value) are accepted,
    "compress": Compressed output under binary format,
    "decompress": Decompressed output (compared with the original dataset),
    "interval": Decompressed interval between two consecutive records,
    "algorithm": {
        "name": Algorithm name,
        "error": Maximum individual error threshold,
        "...": Other algorithm-specific configurations
    }
}
```

- Configuration templates are provided in ``conf/template/``.
- Valid algorithm-specific configurations beyond those predefined in the templates can be identified by inspecting the corresponding code in ``src/C++/*``.


### Execution
To execute an algorithm, run the following command with the corresponding configuration file.

```bash
$ script/run.sh <CONFIG_FILE>
```

### Statistical Output

Statistical output of each execution is appended to the ``experiments.csv`` file, including:

| Column | Description |
|:-------------|:-----------|
| ``Dataset`` | Dataset used during execution. |
| ``Algorithm`` | Algorithm used for compression and decompression. |
| ``Error`` | Maximum allowable individual error threshold. |
| ``Compression ratio`` | Ratio of original data size to compressed data size. |
| ``mse`` | Mean squared error of the reconstructed data. |
| ``rmse`` | Root mean squared error of the reconstructed data. |
| ``mae`` | Mean absolute error of the reconstructed data. |
| ``snr`` | Signal-to-noise ratio of the reconstructed data. |
| ``psnr`` | Peak signal-to-noise ratio of the reconstructed data. |
| ``max_e`` | Maximum error of the reconstructed data. |
| ``min_e`` | Minimum error of the reconstructed data. |
| ``max_vsz`` | Maximum virtual memory size (VSZ) used during execution. |
| ``max_rss`` | Maximum resident set size (RSS) used during execution. |
| ``c_time`` | Average compression time per data point. |
| ``c_avg_latency`` | Average compression latency per data point. |
| ``c_max_latency`` | Maximum compression latency observed. |
| ``d_time`` | Average decompression time per data segment. |
| ``max_d_latency`` | Maximum decompression latency observed. |


## Algorithms
All algorithms from prior work that we have re-implemented are listed below. For the sake of consistency, our naming conventions might look different from those used by the original authors, since some algorithms were not explicitly named in the original papers.

### Piecewise constant approximation:
This family partitions a data stream into multiple segments, with each represented by a **constant value**.
* `pmc` : Capturing Sensor-generated Time Series With Quality Guarantees. (Link: https://ieeexplore.ieee.org/document/1260811)
* `hybrid-pca` : Improved Piecewise Constant Approximation Method for Compressing Data Streams (Link: https://ieeexplore.ieee.org/document/8934460).

### Piecewise linear approximation:
This family partitions a data stream into multiple segments, with each represented by a **linear line**.
* `swing-filter` : Online Piece-wise Linear Approximation of Numerical Streams with Precision Guarantees. (Link: https://dl.acm.org/doi/abs/10.14778/1687627.1687645)
* `slide-filter` : Online Piece-wise Linear Approximation of Numerical Streams with Precision Guarantees. (Link: https://dl.acm.org/doi/abs/10.14778/1687627.1687645) 
* `optimal-pla` : Maximum error-bounded Piecewise Linear Representation for Online Stream Approximation. (Link: https://dl.acm.org/doi/10.1007/s00778-014-0355-0)
* `cov-pla` : Streaming Piecewise Linear Approximation for Efficient Data Management in Edge Computing. (Link: https://dl.acm.org/doi/10.1145/3297280.3297552)
* `conn-I-pla` : An Improved Algorithm for Segmenting Online Time Series with Error Bound Guarantee. (Link: https://link.springer.com/article/10.1007/s13042-014-0310-9)
* `semi-optimal-pla` : An Optimal Online Semi-Connected PLA Algorithm With Maximum Error Bound. (Link: https://ieeexplore.ieee.org/document/9039677)
* `semi-mixed-pla` : An Online PLA Algorithm with Maximum Error Bound for Generating Optimal Mixed‐segments. (Link: https://link.springer.com/article/10.1007/s13042-019-01052-y)
* `mix-piece` : Flexible Grouping of Linear Segments for Highly Accurate Lossy Compression of Time Series Data. (Link: https://link.springer.com/article/10.1007/s00778-024-00862-z)
* `ioriented-pla` : Ours.

### Piecewise polynomial approximation:
This family partitions a data stream into multiple segments, with each represented by a **polynomial of degree k**, where k is a predefined hyperparameter.
* `cached-normal-equation` : Fast Piecewise Polynomial Fitting of Time-Series Data for Streaming Computing. (Link: https://ieeexplore.ieee.org/document/9016024)

### Model selection:
This family partitions a data stream into multiple segments, each represented by **the most suitable model**, which is selected from a set of candidate functions.
* `adaptive-approximation` : An Adaptive Algorithm for Online Time Series Segmentation with Error Bound Guarantee. (Link: https://dl.acm.org/doi/10.1145/2247596.2247620)
* `smart-grid-compression` : A Time-series Compression Technique and Its Application to The Smart Grid. (Link: https://link.springer.com/article/10.1007/s00778-014-0368-8)
* `adapt-ppa` : Ours.

## Contact

For questions, issues, or further information, please contact:

- **Huan** — huan@hcmut.edu.vn  