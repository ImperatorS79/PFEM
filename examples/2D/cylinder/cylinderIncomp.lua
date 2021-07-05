Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 15,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.01,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.55,
		addOnFS = true,
		deleteFlyingNodes = false,
		exclusionZones = {},
		boundingBox = {-0.1, -0.52, 1.05, 0.52},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/2D/cylinder/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.05,
			whatToWrite = {"p", "ke", "velocity"},
			writeAs = "NodesElements" 
		},
		{
			kind = "Mass",
			outputFile = "mass.txt",
			timeBetweenWriting = 1 
		}
	},
	
	Material = {
		mu = 0.015,
		rho = 1000,
		gamma = 0
	},
	
	IC = {
		WallsFixed = true,
		CylinderWallFixed = true,
		FluidInputFixed = true
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.01,
		initialDT = 0.01,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			gammaFS = 0,
			bodyForce = {0, 0},
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0}
end


function Problem.IC:initFluidInputStates(pos)
	return {0, 0, 0}
end

function Problem.Solver.MomContEq.BC:WallsV(pos, t)
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:CylinderWallV(pos, t)
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:FluidInputV(pos, t)
	local v0 = 0.01*4
	DT = 2
	if(t > DT) then
		return{v0*(0.375^4 - pos[2]^4)/0.375^4, 0}
	else
		return {(-2*v0*(t/DT)^3 + 3*v0*(t/DT)^2)*(0.375^4 - pos[2]^4)/0.375^4, 0}
	end
end