# ECA Demo

ECA Demo is a MATLAB repository to share the code and a small-scale example of the simulator and reconstruction engine described in "Enhancement-Constrained Acceleration: a Reconstruction Framework for Breast DCE-MRI" [**ref**].

If you first contact <tyo8teasley@uchicago.edu> or <fdp@uchicago.edu> for datasets, you will also be able to use the code in this repository to reproduce the results in [**ref**].

More detail on the reconstruction algorithm and input/output data strucutres can be found in **Algorithms_and_data_structures.pdf**.

## Installation

Zip and download ECA_demo and extract into a location on the MATLAB path. To keep only the simulation and reconstruction kernels, download the contents of **ECA_demo/Demo/demo_fctns**, ignoring **AP Phantom Functions** and **Paths**.

## Usage

```matlab
>> run_demo(1);  %% 0, 1, 2: demo sampling path options
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)
