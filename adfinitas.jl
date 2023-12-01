using Distributed, DistributedArrays, ProgressMeter;

mutable struct adfinitasVector
    x::BigFloat
    y::BigFloat
    z::BigFloat
    adfinitasVector(x, y, z) = new(x, y, z)
end

mutable struct adfinitasTrack
    name::String
    time::Vector{BigFloat}
    position::Vector{adfinitasVector}
    velocity::Vector{adfinitasVector}
    adfinitasTrack(name) = new(name, BigFloat[], adfinitasVector[], adfinitasVector[])
end

mutable struct adfinitasTrackArray
    name::String
    time::Vector{BigFloat}
    x::Vector{BigFloat}
    y::Vector{BigFloat}
    z::Vector{BigFloat}
    adfinitasTrackArray(name) = new(name, BigFloat[], BigFloat[], BigFloat[], BigFloat[])
end

mutable struct adfinitasObject
    name::String
    mass::BigFloat
    position::adfinitasVector
    velocity::adfinitasVector
    time::BigFloat
    timeStep::BigFloat
    pastPosition::adfinitasVector
    pastVelocity::adfinitasVector
    pastAcceleration::adfinitasVector
    adfinitasObject(name, mass, position, velocity, time) = new(name, mass, position, velocity, time, 0, adfinitasVector(0., 0., 0.), adfinitasVector(0., 0., 0.), adfinitasVector(0., 0., 0.))
end

mutable struct adfinitasParameters
    gravitationalConstant::BigFloat
    timeStep::BigFloat
    angleResolution::BigFloat
    # updateFunction::Function
    adfinitasParameters(gravitationalConstant, timeStep, angleResolution) = new(gravitationalConstant, timeStep, angleResolution)
end

function adfinitasGetArray(vectorObject::adfinitasVector)::Vector{BigFloat}
    x::BigFloat = vectorObject.x;
    y::BigFloat = vectorObject.y;
    z::BigFloat = vectorObject.z;
    return [x, y, z];
end

function adfinitasGetVector(arrayObject::Vector{BigFloat})::adfinitasVector
    return adfinitasVector(arrayObject[1], arrayObject[2], arrayObject[3]);
end

function adfinitasInitObjects(configure::adfinitasParameters, objectsArray::Vector{adfinitasObject})::Vector{adfinitasObject}
    deltaOmega::BigFloat = configure.angleResolution;
    for subObject in objectsArray
        position::Vector{BigFloat} = adfinitasGetArray(subObject.position);
        velocity::Vector{BigFloat} = adfinitasGetArray(subObject.velocity);
        orbitRadius::BigFloat = sqrt(sum(abs.(position) .^ 2));
        if orbitRadius == 0
            subObject.timeStep = configure.timeStep;
        else
            orbitvelocity::BigFloat = sqrt(sum(abs.(velocity) .^ 2));
            subObject.timeStep = deltaOmega * orbitRadius / orbitvelocity;
        end
    end
    return objectsArray;
end

function adfinitasInitSpecifiedObject(configure::adfinitasParameters, object::adfinitasObject, angleCenterObject::adfinitasObject)::adfinitasObject
    initedObj::adfinitasObject = object;
    deltaOmega::BigFloat = configure.angleResolution;

    centerPosition::Vector{BigFloat} = adfinitasGetArray(angleCenterObject.position);
    position::Vector{BigFloat} = adfinitasGetArray(initedObj.position);
    velocity::Vector{BigFloat} = adfinitasGetArray(initedObj.velocity);
    orbitRadius::BigFloat = sqrt(sum(abs.(position .- centerPosition) .^ 2));
    if orbitRadius == 0
        initedObj.timeStep = configure.timeStep;
    else
        orbitvelocity::BigFloat = sqrt(sum(abs.(velocity) .^ 2));
        initedObj.timeStep = deltaOmega * orbitRadius / orbitvelocity;
    end
    
    return initedObj;
end

function adfinitasSubAcceleration(configure::adfinitasParameters, motionObject::adfinitasObject, sourceObject::adfinitasObject)::adfinitasVector

    sourcePosition::Vector{BigFloat} = adfinitasGetArray(sourceObject.position);
    motionPosition::Vector{BigFloat} = adfinitasGetArray(motionObject.position);
    distance::BigFloat = sqrt(sum(abs.(sourcePosition .- motionPosition) .^ 2));

    accelerationMod::BigFloat = sourceObject.mass / (distance ^ 2) * configure.gravitationalConstant;
    accelerationArray::Vector{BigFloat} = accelerationMod / distance .* (sourcePosition .- motionPosition);
    accelerationVector::adfinitasVector = adfinitasGetVector(accelerationArray);

    return accelerationVector;
end

function adfinitasUpdatePosition!(timeStep::BigFloat, motionObject::adfinitasObject, accelerationVector::adfinitasVector)::adfinitasObject
    velocity::Vector{BigFloat} = adfinitasGetArray(motionObject.velocity);
    position::Vector{BigFloat} = adfinitasGetArray(motionObject.position);
    acceleration::Vector{BigFloat} = adfinitasGetArray(accelerationVector);

    # Euler
    # newPosition::Vector{BigFloat} = position .+ (timeStep .* velocity);

    # Velocity Verlet
    newPosition::Vector{BigFloat} = position .+ (timeStep .* velocity) .+ ((timeStep * timeStep / 2) .* acceleration);

    motionObject.position = adfinitasGetVector(newPosition);

    return motionObject;
end

function adfinitasUpdateVelocity!(timeStep::BigFloat, motionObject::adfinitasObject, accelerationVector::adfinitasVector)::adfinitasObject
    velocity::Vector{BigFloat} = adfinitasGetArray(motionObject.velocity);
    # position::Vector{BigFloat} = adfinitasGetArray(motionObject.position);
    acceleration::Vector{BigFloat} = adfinitasGetArray(accelerationVector);

    # Euler
    # newVelocity::Vector{BigFloat} = velocity .+ (timeStep .* acceleration);

    # Velocity Verlet
    lastAcceleration::Vector{BigFloat} = adfinitasGetArray(motionObject.pastAcceleration);
    newVelocity::Vector{BigFloat} = velocity .+ ((timeStep / 2) .* (acceleration .+ lastAcceleration));

    motionObject.velocity = adfinitasGetVector(newVelocity);
    motionObject.pastAcceleration = accelerationVector;

    return motionObject;
end

function adfinitasAcceleration(configure::adfinitasParameters, motionObject::adfinitasObject, objectsArray::Vector{adfinitasObject})::adfinitasVector
    accelerationArray::Vector{BigFloat} = zeros(BigFloat, 3);
    for sourceObject in objectsArray
        if motionObject.name != sourceObject.name
            accelerationArray += adfinitasGetArray(adfinitasSubAcceleration(configure, motionObject, sourceObject));
        end
    end
    accelerationVector::adfinitasVector = adfinitasGetVector(accelerationArray);
    return accelerationVector;
end

function adfinitasOneStep(configure::adfinitasParameters, objectsArray::Vector{adfinitasObject}, isPosition::Bool)::Vector{adfinitasObject}

    distributedClonedObjectsArray::DArray = distribute([adfinitasObject[] for _ in procs()]);

    @sync @distributed for motionObjectLoop in objectsArray
        motionObject::adfinitasObject = deepcopy(motionObjectLoop);
        sysTimeStep::BigFloat = configure.timeStep;
        objTimeStep::BigFloat = motionObject.timeStep;
        aimTime::BigFloat = motionObject.time + sysTimeStep;
        timeStep::BigFloat = 0.;
        loopTimes::Int128 = 0;
        if sysTimeStep > objTimeStep
            loopTimes = floor(sysTimeStep / objTimeStep) + 1;
            timeStep = objTimeStep;
        else
            timeStep = sysTimeStep;
            loopTimes = 1;
        end
        for loopIndex in 1 : loopTimes
            if loopIndex == loopTimes
                finishTime::BigFloat = motionObject.time + timeStep;
                if finishTime > aimTime
                    timeStep = aimTime - motionObject.time;
                end
            end
            accelerationVector::adfinitasVector = adfinitasVector(0., 0., 0.);
            if isPosition
                if sum(abs.(adfinitasGetArray(motionObject.pastAcceleration))) == 0.
                    accelerationVector = adfinitasAcceleration(configure, motionObject, objectsArray);
                else
                    accelerationVector = motionObject.pastAcceleration;
                end
                motionObject = adfinitasUpdatePosition!(timeStep, motionObject, accelerationVector);
                motionObject.time += timeStep;
            else
                accelerationVector = adfinitasAcceleration(configure, motionObject, objectsArray);
                motionObject = adfinitasUpdateVelocity!(timeStep, motionObject, accelerationVector);
            end
        end
        push!(localpart(distributedClonedObjectsArray)[1], motionObject);
    end

    objectsArrayFinish::Vector{adfinitasObject} = [];
    for singleObject in distributedClonedObjectsArray
        append!(objectsArrayFinish, singleObject);
    end

    distributedClonedObjectsArray = distribute([nothing for _ in procs()]);

    return objectsArrayFinish;
end

function adfinitasGetTrack(objectName::String, objectsArray::Vector{adfinitasObject})::adfinitasTrack
    track::adfinitasTrack = adfinitasTrack(objectName);
    for subObj in objectsArray
        if subObj.name == objectName
            time::BigFloat = subObj.time;
            position::adfinitasVector = subObj.position;
            velocity::adfinitasVector = subObj.velocity;
            append!(track.time, time);
            append!(track.position, [position]);
            append!(track.velocity, [velocity]);
        end
    end
    return track;
end

function adfinitasConvertTrackArray(track::adfinitasTrack)::Tuple{adfinitasTrackArray, adfinitasTrackArray}
    trackPositionArray::adfinitasTrackArray = adfinitasTrackArray(track.name);
    trackVelocityArray::adfinitasTrackArray = adfinitasTrackArray(track.name);
    index::Int128 = 1;
    for time in track.time
        append!(trackPositionArray.time, time);
        append!(trackPositionArray.x, track.position[index].x);
        append!(trackPositionArray.y, track.position[index].y);
        append!(trackPositionArray.z, track.position[index].z);
        append!(trackVelocityArray.time, time);
        append!(trackVelocityArray.x, track.velocity[index].x);
        append!(trackVelocityArray.y, track.velocity[index].y);
        append!(trackVelocityArray.z, track.velocity[index].z);
        index += 1;
    end
    return trackPositionArray, trackVelocityArray;
end

function adfinitasRun(configure::adfinitasParameters, simulationTime::BigFloat, objsArray::Vector{adfinitasObject})::Vector{adfinitasObject}
    objsMotionList::Vector{adfinitasObject} = [];
    simulationLength::Int128 = convert(Int128, floor(simulationTime / configure.timeStep));

    @showprogress for times in 1 : simulationLength
    # for times in 1 : simulationLength
        objsArray = adfinitasOneStep(configure, objsArray, true);
        objsArray = adfinitasOneStep(configure, objsArray, false);
        append!(objsMotionList, objsArray);
    end

    return objsMotionList;
end

function adfinitasEmptyObjectArray()::Vector{adfinitasObject}
    emptyObjArray::Vector{adfinitasObject} = [];
    return emptyObjArray;
end

function adfinitasParallelStart(nProcs::Int64)
    currentProcs::Int128 = nprocs();
    if currentProcs < nProcs
        addprocs(nProcs - currentProcs);
    end
end

function adfinitasParallelStop()
    for proc in procs()
        rmprocs(proc);
    end
end

function adfinitasDistributedStart(server::String, procs_n::Int64, shell::Symbol, workingDir::String)
    # server(SSH): "user@address:port";
    # procs_n: remote procs number;
    # shell(Symbol): :csh, :bash, :zsh...
    addprocs([(server, procs_n)], tunnel=true, shell=shell, dir=workingDir, exename="julia");
end
