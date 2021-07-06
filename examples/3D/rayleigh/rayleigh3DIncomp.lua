Problem = {
    id = "Boussinesq",
	simulationTime = 50,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.1,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.35,
		addOnFS = true,
		deleteFlyingNodes = false,
		boundingBox = {-0.05, -0.05, -0.05, 4.05, 4.05, 100},
		exclusionZones = {},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/3D/rayleigh/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.5,
			whatToWrite = {"T", "ke", "p"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1e-3,
		gamma = 0,
		rho = 1000,
		k = 0.6,
		cv = 1,
		alpha = 69e-6,
		Tr = 650,
		DgammaDT = 0,
		h = 0,
		Tinf = 300,
		epsRad = 0
	},
	
	IC = {
		LateralFixed = true,
		TopFixed = true,
		BottomFixed = true
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.005,
		maxRemeshDT = -1,
		initialDT = 0.005,
		solveHeatFirst = true,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, 0, -9.81},
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
			
			}
		}
	}
}

function Problem.IC:initStates(pos)
    local g = -Problem.Solver.MomContEq.bodyForce[3]
	local rho = Problem.Material.rho
	local p = rho*g*(1 - pos[3])

	return {0, 0, 0, p, 650}
end

function Problem.IC:initTopStates(pos)
	local g = -Problem.Solver.MomContEq.bodyForce[3]
	local rho = Problem.Material.rho
	local p = rho*g*(1 - pos[3])

	return {0, 0, 0, p, 300}
end

function Problem.IC:initBottomStates(pos)
	local g = -Problem.Solver.MomContEq.bodyForce[3]
	local rho = Problem.Material.rho
	local p = rho*g*(1 - pos[3])

	return {0, 0, 0, p, 1000}
end

function Problem.Solver.MomContEq.BC:LateralV(pos, t)
	return {0, 0, 0}
end

function Problem.Solver.MomContEq.BC:TopV(pos, t)
	return {0, 0, 0}
end

function Problem.Solver.MomContEq.BC:BottomV(pos, t)
	return {0, 0, 0}
end

function Problem.Solver.HeatEq.BC:TopT(pos, t)
	return {300}
end

function Problem.Solver.HeatEq.BC:BottomT(pos, t)
	return {1000}
end
