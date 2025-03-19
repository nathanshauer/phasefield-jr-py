# phasefield-jr-py
#### A simple one file Python project to run 2D phase field problems with linear quadrilaterals

The phase field method is a powerful tool for fracture analysis. However, it introduces certain challenges that are not encountered in traditional finite element analysis. With this in mind, this code was developed for educational purposes, providing a single-file implementation to help researchers familiarize themselves with the fundamentals of phase field analysis. It also serves as a reference for verifying their own code.

Check out the sister code in C++ on:
[https://github.com/nathanshauer/phasefield-jr](https://github.com/nathanshauer/phasefield-jr)

For more information about me or to get in touch, please visit my website:
[www.nathanshauer.com](http://www.nathanshauer.com)

## Configuration

The code has been tested on macOS and Ubuntu. 

### Installing Python, numpy, scipy and matplotlib

To run the code, you need to have Python, numpy and matplotlib installed on your system. Scipy is also needed if using the examples that adopt sparse matrices (`example1.py` and `example2.py`). Follow these steps to install them:

1. **Install Python**: If you don't have Python installed, download and install it from the [official website](https://www.python.org/downloads/). You can also download python using package manager in Linux or macports/homebrew in macOS. **The code was tested using Python 3.12**

2. **Install numpy and matplotlib**: Open a terminal and run the following commands to install numpy and matplotlib using pip:

```sh
pip install numpy scipy matplotlib
```

3. **Verify installation**: To ensure that the packages are installed correctly, you can run the following commands in a Python shell:

```python
import numpy
import scipy
import matplotlib
print(numpy.__version__)
print(scipy.__version__)
print(matplotlib.__version__)
```

If the versions are printed without any errors, the installation was successful.

## Running the code

This repository includes two sets of examples: one using sparse matrices and another using dense matrices for the global stiffness matrix. 

1. **Sparse Matrix Examples**: 
  - `example1_sparse.py`: Simulates a bar under tension using sparse matrices for improved scalability.
  - `example2_sparse.py`: Simulates a single-edge notch plate under tension using sparse matrices.

  Sparse matrices are more memory-efficient and computationally scalable, especially for large problems. However, they introduce additional complexities in data structure management.

2. **Dense Matrix Examples**:
  - `example1_dense.py`: A dense matrix version of the bar under tension simulation.
  - `example2_dense.py`: A dense matrix version of the single-edge notch plate simulation.

  Dense matrices are easier to understand and implement, making these examples more suitable for educational purposes or smaller problems.

Both versions produce comparable results.

To run any example, open a terminal, navigate to the project directory, and execute the following command:

```sh
python <example_file>.py
```

Replace `<example_file>` with the desired example script (e.g., `example1_sparse.py` or `example2_dense.py`).

Both examples 1 and 2 will generate output files in the `output` directory, which can be visualized using ParaView as described in the "Output in Paraview using vtk files" section.

## Output in Paraview using vtk files
To visualize the output in ParaView using VTK files, follow these steps:

1. **Generate VTK files**: After running the examples, VTK files will be generated in the `output` directory. These files contain the simulation results and can be visualized using ParaView.

2. **Install ParaView**: If you don't have ParaView installed, you can download it from the [official website](https://www.paraview.org/download/).

3. **Open ParaView**: Launch ParaView on your system.

4. **Load VTK files**:
  - Click on `File` > `Open`.
  - Navigate to the `output` directory of the project.
  - Select the VTK file you want to visualize (e.g., `output_ex1_#.vtk` or `output_ex2_#.vtk`).
  - Click `Apply` to load the data.

5. **Visualize the data**: Use the various visualization tools in ParaView to explore the simulation results. You can adjust the display properties, apply filters, and create animations to better understand the phase field analysis.

## Quantitative analyses 

The project includes quantitative analysis of the simulation results directly within the Python scripts. At the end of each run, the analysis is performed, and the results are saved in the `output` directory.

1. **Analysis of Bar under Tension (example1)**
  - This analysis is performed at the end of the `example1` simulation.
  - The script generates a plot of `sigma/sigma_peak x time`, where `sigma_peak` is calculated analytically.
  - The plot is saved as `example1_timexsigma.png` in the `output` directory.

2. **Analysis of Single-edge Notch Plate under Tension (example2)**
  - This analysis is performed at the end of the `example2` simulation.
  - The script generates a plot of `Reaction force x imposed displacement`.
  - The plot is saved as `example2_pdelta.png` in the `output` directory.

To view the analysis results, navigate to the `output` directory after running the simulations and open the respective PNG files.