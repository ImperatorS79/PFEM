Problem = {
    id = "Boussinesq",
	simulationTime = 50,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.02,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		addOnFS = true,
		deleteFlyingNodes = false,
		boundingBox = {-0.05, -0.05, 1.05, 1.05},
		exclusionZones = {},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/2D/thermalConv/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 1,
			whatToWrite = {"T", "ke", "p"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1e-3,
		rho = 1000,
		k = 0.6,
		alpha = 69e-6,
		Tr = 300,
		cv = 4.186,
		gamma = 0,
		DgammaDT = 0,
		h = 0,
		Tinf = 300,
		epsRad = 0
	},
	
	IC = {
		TopFixed = true,
		BottomFixed = true,
		LeftFixed = true,
		RightFixed = true
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.2,
		initialDT = 0.2,
		solveHeatFirst = true,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, -9.81},
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
	return {0, 0, 0, 300}
end

function Problem.IC:initLeftStates(pos)
	return {0, 0, 0, 310}
end

function Problem.IC:initRightStates(pos)
	return {0, 0, 0, 290}
end

function Problem.Solver.HeatEq.BC:TopQ(pos, t) 
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:TopV(pos, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:BottomQ(pos, t) 
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:BottomV(pos, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:LeftT(pos, t) 
	return {310}
end

function Problem.Solver.MomContEq.BC:LeftV(pos, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:RightT(pos, t) 
	return {290}
end

function Problem.Solver.MomContEq.BC:RightV(pos, t) 
	return {0, 0}
end
