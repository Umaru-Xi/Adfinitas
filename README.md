# Adfīnitās

Julia N-body simulation package, support parallel and distributed computing.  

* * *  

Package dependency: Distributed, DistributedArrays, ProgressMeter;

Example dependenchy: Jupyter Notebook, PyPlot, PyCall;

* * *  

Example: examples/examples.ipynb

1. Solar system up to Jupiter;  
    a. Plot orbits on 3D;  
        ![solar system orbits 3D](examples/Result/position.svg.png)  
    b. Animate planets and solar's tracks;  
        ![animate solar system orbits](examples/Result/animePosition.gif)  
    c. Plot the radial velocity of Sun and Mercury;  
        ![radial velosity of Sun](examples/Result/SunRadialVelocity.svg.png)  
        ![radial velosity of Mercury](examples/Result/MercuryRadialVelocity.svg.png)  
2. Sun-Earth-Moon system, but Moon has velocity on z-axis;  
    a. 3D orbits plot;  
        ![Sun-Earth-Moon system orbits 3D](examples/Result/MoonPosition.svg.png)  
3. N-Body simulation with a galaxy;  
    a. 100-body animate on 2D;  
        ![100-body animation 2D](examples/Result/anime2DGalaxyPosition.gif)  
    b. 1000-body animate on 3D;  
        ![1000-body animation 3D](examples/Result/animeGalaxyPosition.gif)  
