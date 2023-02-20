# Getting Started
## Installation and Setup
This should be a pretty easy install. These instructions assume you have MATLAB installed. The scripts were written in R2022a, but other versions should work too.
1. Download these scripts (into the same folder).
2. Create a data folder with the `.edf` files in it. This data folder could be anywhere, but it's convinient to put it in the same folder as the scripts.
3. Install the Matlab version of (chronux)[http://chronux.org/]. In this case, this means downloading the library so there's a folder somewhere with the chronux code.
4. Modify `get_directory_info.m` to point to the appropriate directories on your machine.
5. Modify `get_clip_metadata.m` to reflect the clips in your data directory (this file will change based on your experiment).

The data folder should have this structure:
```
[data folder]/
    baseline/
        [animal]_[temp].edf
        ...
    seizure/
        [animal]_[temp].edf
        ...
```
So for example, the setup the scripts are currently configured for is:
```
data_folder/ (the data folder)
    baseline/
        B_37.edf
        C_37.edf
        D_37.edf
        E_37.edf
        F_37.edf
        L_37.edf
        L_42.edf
        N_37.edf
        N_42.edf
        O_37.edf
    seizure/
        F_37.edf
        O_42.edf
generated/ (the output folder)
calculate_comodulogram.m
...
single_channel.m

../../../Documents/MATLAB/chronux_2_11 (the chronux folder)
```

## Usage
The two main scripts you'll want to run will be `single_channel.m` and `multi_channel.m`. Everything else (except `comodulogram_test.m`) supports those two files. `single_channel.m` considers the EEG's individually, while `multi_channel.m` considers groups of them. 
Conceptually, the single channel analysis is easier to understand, but I spent more time in `multi_channel.m`, so that's better thought-out and documented. 
Many of the cells in the two scripts paralell each other, so if something doesn't make sense in one, check the other.

# File explanations
Note that get_ functions work with the filesystem while calculate_ functions are pure functions. All other files are scripts.
* calculate_comodulogram: a function that to create comodulograms given a signal and frequency bands
* calculate_comodulogram_stack: a function similar to calculate_comodulogram, but it first breaks the signal into chunks and returns the comodulogram for each chunck
* comodulogram_test: a script that walks through the process of developing a comodulogram and demonstrates some of their limitations
* get_clip_metadata: a function that holds the information about the clips that get analysed in the _channel scripts.
* get_directory_info: a function that stores the paths to important files; it's basically a configuration file
* get_lfp: this file reads and minimally processes `.edf` files into LFP signals
* multi_clip.m: this script analyses a set of clips in paralell (formerly multi_channel.m)
* single_clip.m: this script analyses clips one at a time (formerly single_channel.m)

# Contact
The original author of these scripts is Jonathan Gould (jngould@umich.edu).

# Acknoweledgements
* Noor Daddo
* Julie Ziobro
* Jack Parent
* Wei Nu 
* Dan Jacklick
* The University of Michigan Neural Graduate Program

# License
Expat License

Copyright (c) 2023 Jonathan Gould

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
