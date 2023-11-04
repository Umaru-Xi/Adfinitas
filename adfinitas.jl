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
    adfinitasObject(name, mass, position, velocity, time) = new(name, mass, position, velocity, time)
end

mutable struct adfinitasParameters
    gravationalConstant::BigFloat
    timeStep::BigFloat
    angleResolution::BigFloat
    adfinitasParameters(gravationalConstant, timeStep, angleResolution) = new(gravationalConstant, timeStep, angleResolution)
end

function adfinitasGetArray(vectorObject::adfinitasVector)
    x = vectorObject.x;
    y = vectorObject.y;
    z = vectorObject.z;
    return [x, y, z];
end

function adfinitasGetVector(arrayObject::Vector{BigFloat})
    return adfinitasVector(arrayObject[1], arrayObject[2], arrayObject[3]);
end

function adfinitasInitObjects(configure::adfinitasParameters, objectsArray::Vector{Any})
    deltaOmega = configure.angleResolution;
    for subObject in objectsArray
        position = adfinitasGetArray(subObject.position);
        velocity = adfinitasGetArray(subObject.velocity);
        orbitRadius = sqrt(sum(abs.(position) .^ 2));
        if orbitRadius == 0
            subObject.timeStep = configure.timeStep;
        else
            orbitVelosity = sqrt(sum(abs.(velocity) .^ 2));
            subObject.timeStep = deltaOmega * orbitRadius / orbitVelosity;
        end
    end
    return objectsArray;
end

function adfinitasInitSpecifiedObject(configure::adfinitasParameters, object::adfinitasObject, angleCenterObject::adfinitasObject)
    initedObj = deepcopy(object);
    deltaOmega = configure.angleResolution;

    centerPosition = adfinitasGetArray(angleCenterObject.position);
    position = adfinitasGetArray(initedObj.position);
    velocity = adfinitasGetArray(initedObj.velocity);
    orbitRadius = sqrt(sum(abs.(position .- centerPosition) .^ 2));
    if orbitRadius == 0
        initedObj.timeStep = configure.timeStep;
    else
        orbitVelosity = sqrt(sum(abs.(velocity) .^ 2));
        initedObj.timeStep = deltaOmega * orbitRadius / orbitVelosity;
    end
    
    return initedObj;
end

function adfinitasSubAcceleration(configure::adfinitasParameters, motionObject::adfinitasObject, sourceObject::adfinitasObject)
    gravationalConstant = configure.gravationalConstant;

    sourcePosition = adfinitasGetArray(sourceObject.position);
    motionPosition = adfinitasGetArray(motionObject.position);
    distance = sqrt(sum(abs.(sourcePosition .- motionPosition) .^ 2));

    accelerationMod = sourceObject.mass / (distance ^ 2) * gravationalConstant;
    accelerationArray = accelerationMod / distance .* (sourcePosition .- motionPosition);
    accelerationVector = adfinitasGetVector(accelerationArray);

    return accelerationVector;
end

function adfinitasUpdate(timeStep, motionObject::adfinitasObject, accelerationVector::adfinitasVector)
    velocity = adfinitasGetArray(motionObject.velocity);
    position = adfinitasGetArray(motionObject.position);
    acceleration = adfinitasGetArray(accelerationVector);

    newVelocity = velocity .+ (timeStep .* acceleration);
    newPosition = position .+ ((timeStep / 2) .* velocity) .+ ((timeStep / 2) .* newVelocity);

    motionObject.velocity = adfinitasGetVector(newVelocity);
    motionObject.position = adfinitasGetVector(newPosition);

    return motionObject;
end

function adfinitasAcceleration(configure::adfinitasParameters, motionObject::adfinitasObject, objectsArray::Vector{Any})
    accelerationArray = zeros(BigFloat, 3);
    for sourceObject in objectsArray
        if motionObject.name != sourceObject.name
            accelerationArray += adfinitasGetArray(adfinitasSubAcceleration(configure, motionObject, sourceObject));
        end
    end
    accelerationVector = adfinitasGetVector(accelerationArray);
    return accelerationVector;
end

function adfinitasOneStep(configure::adfinitasParameters, objectsArray::Vector{Any})

    clonedObjectsArray = deepcopy(objectsArray);
    distributedClonedObjectsArray = distribute([adfinitasObject[] for _ in procs()]);

    @sync @distributed for motionObject in clonedObjectsArray
        sysTimeStep = configure.timeStep;
        objTimeStep = motionObject.timeStep;
        aimTime = motionObject.time + sysTimeStep;
        if sysTimeStep > objTimeStep
            loopTimes = floor(sysTimeStep / objTimeStep) + 1;
            timeStep = objTimeStep;
        else
            timeStep = sysTimeStep;
            loopTimes = 1;
        end
        for loopIndex in 1 : loopTimes
            if loopIndex == loopTimes
                finishTime = motionObject.time + timeStep;
                if finishTime > aimTime
                    timeStep = aimTime - motionObject.time;
                end
            end
            accelerationVector = adfinitasAcceleration(configure, motionObject, objectsArray);
            motionObject = adfinitasUpdate(timeStep, motionObject, accelerationVector);
            motionObject.time += timeStep;
        end
        push!(localpart(distributedClonedObjectsArray)[1], motionObject);
    end

    objectsArray = [];
    for singleObject in distributedClonedObjectsArray
        append!(objectsArray, singleObject);
    end
    return objectsArray;
end

function adfinitasGetTrack(objectName::String, objectsArray::Vector{Any})
    track = adfinitasTrack(objectName);
    for subObj in objectsArray
        if subObj.name == objectName
            time = subObj.time;
            position = subObj.position;
            velocity = subObj.velocity;
            append!(track.time, time);
            append!(track.position, [position]);
            append!(track.velocity, [velocity]);
        end
    end
    return track;
end

function adfinitasConvertTrackArray(track::adfinitasTrack)
    trackPositionArray = adfinitasTrackArray(track.name);
    trackVelocityArray = adfinitasTrackArray(track.name);
    index = 1;
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

function adfinitasParallelStart(nProcs::Int64)
    currentProcs = nprocs();
    if currentProcs < nProcs
        addprocs(nProcs - currentProcs);
    end
end

function adfinitasParallelStop()
    for proc in procs()
        rmprocs(proc);
    end
end

function adfinitasRun(configure::adfinitasParameters, simulationTime::BigFloat, objsArray::Vector{Any})
    objsMotionList = [];
    simulationLength = convert(Int64, floor(simulationTime / configure.timeStep));

    @showprogress for times in 1 : simulationLength
        objsArray = adfinitasOneStep(configure, objsArray);
        append!(objsMotionList, objsArray);
    end

    return objsMotionList;
end