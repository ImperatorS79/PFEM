Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 4,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.0073,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-2, -1, 12, 100},
		addOnFS = true,
		deleteFlyingNodes = true,
		exclusionZones = {},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/2D/dropFall/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.01,
			whatToWrite = {"p", "ke"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1e-3,
		gamma = 0,
		K0 = 2200000,
		K0p = 7.6,
		rhoStar = 1000
	},
	
	IC = {
		BoundaryFixed = true,
	},
	
	Solver = {
	    id = "CDS_dpdt",
		adaptDT = true,
		maxDT = 0.0001,
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
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0, Problem.Material.rhoStar, 0, 0}
end

function Problem.Solver.MomEq.BC:BoundaryV(pos, t)
	return {0, 0}
end
