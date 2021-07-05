Problem = {
    id = "Boussinesq",
	simulationTime = 50,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.02,
		alpha = 1.25,
		omega = 0.7,
		gamma = 0.6,
		addOnFS = true,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		boundingBox = {-0.05, -0.05, 4.05, 100},
		exclusionZones = {},
		mshFile = "examples/2D/rayleigh/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.1,
			whatToWrite = {"T", "ke", "p"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1e-3,
		rho = 1000,
		k = 0.6,
		alpha = 69e-6,
		Tr = 650,
		cv = 1,
		gamma = 0,
		DgammaDT = 0,
		h = 5,
		Tinf = 300,
		epsRad = 0
	},
	
	IC = {
		BottomFixed = true,
		LeftFixed = true,
		RightFixed = true
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.01,
		initialDT = 0.01,
		solveHeatFirst = true,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			computePres = true,
			bodyForce = {0, -9.81},
			gammaFS = 1,
			residual = "Ax_f",
			BC = {

			}
		},
		
		HeatEq = {
			minRes = 1e-6,
			maxIter = 10,
			residual = "Ax_f",
			BC = {
				FreeSurfaceQh = true
			}
		}
	}
}

function Problem.IC:initStates(pos)
	local g = -Problem.Solver.MomContEq.bodyForce[2]
	local rho = Problem.Material.rho
	local p = rho*g*(1 - pos[2])

	return {0, 0, p, 650}
end

function Problem.IC:initTopStates(pos)
	local g = -Problem.Solver.MomContEq.bodyForce[2]
	local rho = Problem.Material.rho
	local p = rho*g*(1 - pos[2])

	return {0, 0, p, 300}
end

function Problem.IC:initBottomStates(pos)
	local g = -Problem.Solver.MomContEq.bodyForce[2]
	local rho = Problem.Material.rho
	local p = rho*g*(1 - pos[2])
	
	return {0, 0, p, 1000}
end

function Problem.Solver.HeatEq.BC:TopT(pos, initPos, states, t) 
	return {300}
end

function Problem.Solver.MomContEq.BC:TopV(pos, initPos, states, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:BottomT(pos, t) 
	return {1000}
end

function Problem.Solver.MomContEq.BC:BottomV(pos, t) 
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:LeftV(pos, t) 
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:RightV(pos, t) 
	return {0, 0}
end
