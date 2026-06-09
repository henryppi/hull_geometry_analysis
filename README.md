// hull_geometry_analysis
// analyze step/stl hull shapes

// install
conda create --name conda_env_para python=3.12
conda activate conda_env_para
conda install -c conda-forge scipy numpy matplotlib vtk svgpathtools
conda install -c conda-forge pythonocc-core
conda install pip
pip install bezier

// fix bezier lib
pip install bezier # conda doesn't have bezier
// fix/remove executable stack flag
sudo apt install patchelf
patchelf --clear-execstack ~/sources/thirdparty/miniconda3/envs/pyoccenv/lib/python3.12/site-packages/bezier.libs/libbezier-e5395c70.so.2024.6.20 



// doesnt work 
 conda env export > environment.yml -n conda

// extract images
pdfimages -all 16229.pdf rudder_ 

//extract svg  
pdftocairo -svg 16229.pdf output.svg
