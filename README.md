# 2L-DRM_flood
2L-DRM fortran codes for Kim and Domiguez (2023) flood study
Have modified the [original 2L-DRM](https://github.com/huancui/DRM_2LDRM). Now reading binary formatted input files. For detailed information, refer to the following reference.

Reference: [Dominguez et al. (2020)](https://journals.ametsoc.org/view/journals/hydr/21/1/jhm-d-19-0101.1.xml)

### Input files:
1.	Forcing data:

|  | File names | Description |
| ------ | ------ | ------ |
| A | ET_yyyy_lamb.dat | evapotranspiration (ET) with daily resolution |
| B | P2_yyyy_lamb.dat | precipitation (PP) with daily resolution* |
| C | QWU_yyyy_lamb.dat | precipitable water for the upper slab (upper-slab integrated vapor) with daily resolution | 
| D | QWD_yyyy_lamb.dat | precipitable water for the lower slab (lower-slab integrated vapor) with daily resolution | 
| E | U3U_yyyy_lamb.dat | upper-slab vapor weighted U with 3-hourly resolution |
| F | U3D_yyyy_lamb.dat | lower-slab vapor weighted U with 3-hourly resolution |
| G | V3U_yyyy_lamb.dat | upper-slab vapor weighted V with 3-hourly resolution |
| H | V3D_yyyy_lamb.dat | lower-slab vapor weighted V with 3-hourly resolution |
| I | W3U_yyyy_lamb.dat | upward vertical wind W at 700mb level with 3-hourly resolution |
| J | W3D_yyyy_lamb.dat | downward vertical wind W at 700mb level with 3-hourly resolution |
| K | q700_yyyy_lamb.dat | humidity at 700mb level with 3-hourly resolution |

(*Precipitation is used to weight regional contribution, is not necessary if replaced by PW to weight)

2.	Region setup files:

|  | File names | description |
| ------ | ------ | ------ |
| A | MW_il.dat | list of i index for target region grids in the domain |
| B | MW_jl.dat | list of j index for target region grids in the domain |
| C | ERA5_lcgd.dat | regional masks for each sub-region defined in the domain (specify each sub-region with a different value) |
| D | 700_masklc.dat | masks to indicate if the grid has only the upper slab (elevation is above 700mb) or both slabs (0 – only upper slab, 1 – both slabs) | 


### Output files:

|  | File names | Description |
| ------ | ------ | ------ |
| A* | yyyy_DAYRA_##.txt | ratio contributions from each sub-region to target region averaged from upper and lower slabs | 
| B* | yyyy_DAYRU_##.txt | ratio contributions from each sub-region to target region only for the upper slab |
| C* | yyyy_DAYRL_##.txt | ratio contributions from each sub-region to target region only for the lower slab |
| D | yyyy_trj_U_##.txt | recording information along each back-trajectory for the upper slab |
| E | yyyy_trj_D_##.txt | recording information along each back-trajectory for the lower slab |
| F** | yyyy_8MLCpc_ddd_UPPER2L.dat | ratio contributions from each grid to target region averaged from upper slab | 
| G** | yyyy_8MLCpc_ddd_LOWER2L.dat | ratio contributions from each grid to target region averaged from lower slab | 
| H** | yyyy_8MLCpc_ddd_AVERA2L.dat | ratio contributions from each grid to target region averaged from upper and lower slabs | 

(Regional contributions are weighted by precipitable water)
("##" indicates abbreviation of target region name)
("*" indicates outputs of the subregion simulation)
("**" indicates outputs of the source simulation)
