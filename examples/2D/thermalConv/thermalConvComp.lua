Problem = {
    id = "BoussinesqWC",
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
		k = 0.6,
		alpha = 69e-6,
		Tr = 300,
		cv = 4.186,
		gamma = 0,
		K0 = 2200000,
		K0p = 7.6,
		rhoStar = 1000,
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
	    id = "CDS_dpdt",
		adaptDT = true,
		maxDT = 0.02,
		initialDT = 1e-8,
		securityCoeff = 0.1,
		
		MomEq = {
			bodyForce = {0, -9.81},
			BC = {
			
			}
		},
		
		ContEq = {
			stabilization = "Meduri",
			BC = {

			}
		},
		
		HeatEq = {
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local rhoStar = Problem.Material.rhoStar
	local g = -Problem.Solver.MomEq.bodyForce[2]

	local rho, p
	rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(1 - pos[2]) + 1)^(1/(K0p - 1))
	p = K0/K0p*((rho/rhoStar)^K0p - 1)

	return {0, 0, p, rho, 0, 0, 300}
end

function Problem.IC:initLeftStates(pos)
    local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local rhoStar = Problem.Material.rhoStar
	local g = -Problem.Solver.MomEq.bodyForce[2]

	local rho, p
	local rho, p
	rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(1 - pos[2]) + 1)^(1/(K0p - 1))
	p = K0/K0p*((rho/rhoStar)^K0p - 1)
	
	return {0, 0, p, rho, 0, 0, 310}
end

function Problem.IC:initRightStates(pos)
	local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local rhoStar = Problem.Material.rhoStar
	local g = -Problem.Solver.MomEq.bodyForce[2]

	local rho, p
	rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(1 - pos[2]) + 1)^(1/(K0p - 1))
	p = K0/K0p*((rho/rhoStar)^K0p - 1)
	
	return {0, 0, p, rho, 0, 0, 290}
end

function Problem.Solver.HeatEq.BC:TopQ(pos, t) 
	return {0, 0}
end

function Problem.Solver.MomEq.BC:TopV(pos, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:BottomQ(pos, t) 
	return {0, 0}
end

function Problem.Solver.MomEq.BC:BottomV(pos, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:LeftT(pos, t) 
	return {310}
end

function Problem.Solver.MomEq.BC:LeftV(pos, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:RightT(pos, t) 
	return {290}
end

function Problem.Solver.MomEq.BC:RightV(pos, t) 
	return {0, 0}
end

