T. Easley, Z. Ren, B. Kim, G. Karczmar, R.F. Barber, F.D. Pineda.

Nov. 2, 2020

# ECA Demo

ECA Demo is a MATLAB repository to share the code and a small-scale example of the simulator and reconstruction engine described in "Enhancement-Constrained Acceleration: a Reconstruction Framework for Breast DCE-MRI" [**ref**].

The code in this repository may be run on the data found in "Case Data" to reproduce the results in [**ref**].

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
