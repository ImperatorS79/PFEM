Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 6,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.04,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		addOnFS = true,
		deleteFlyingNodes = false,
		boundingBox = {-0.01, -0.01, 5, 1.01},
		exclusionZones = {},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/2D/pipe/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.5,
			whatToWrite = {"p", "u"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 200,
		gamma = 0,
		K0 = 2200000,
		K0p = 7.6,
		rhoStar = 1000
	},
	
	IC = {
		BoundaryFixed = true,
		FluidInputFixed = true,
	},
	
	Solver = {
	    id = "CDS_dpdt",
		adaptDT = true,
		maxDT = 0.005,
		initialDT = 1e-8,
		securityCoeff = 0.1,
		
		MomEq = {
			bodyForce = {0, 0},
			BC = {
			
			}
		},
		
		ContEq = {
			stabilization = "Meduri",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
    local rhoStar = Problem.Material.rhoStar
	return {1, 0, 0, rhoStar, 0, 0}
end

function Problem.IC:initBoundaryStates(pos)
	local rhoStar = Problem.Material.rhoStar
	return {0, 0, 0, rhoStar, 0, 0}
end

function Problem.Solver.MomEq.BC:BoundaryV(pos, t)
	return {0, 0}
end

function Problem.Solver.MomEq.BC:FluidInputV(pos, t)
	return {0, 0}
end