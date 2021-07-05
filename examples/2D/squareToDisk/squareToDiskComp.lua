Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 4,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.05,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.4,
		boundingBox = {-2, -2, 2, 2},
		exclusionZones = {},
		addOnFS = true,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/2D/squareToDisk/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.1,
			whatToWrite = {"p", "ke"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1.0,
		gamma = 1.9,
		K0 = 220000,
		K0p = 7.6,
		rhoStar = 100
	},
	
	IC = {
		BoundaryFixed = true,
	},
	
	Solver = {
	    id = "CDS_dpdt",
		adaptDT = true,
		maxDT = 0.001,
		initialDT = 1e-8,
		securityCoeff = 0.025,
		
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
	local K0p = Problem.Material.K0p
	local K0 = Problem.Material.K0
	
	return {0, 0, 0, rhoStar, 0, 0}
end
