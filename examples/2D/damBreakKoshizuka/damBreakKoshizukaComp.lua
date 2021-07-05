Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 4,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.0073,
		alpha = 1.2,
		omega = 0.6,
		gamma = 0.6,
		addOnFS = true,
		deleteFlyingNodes = false,
		boundingBox = {-0.01, -0.01, 0.594, 100},
		laplacianSmoothingBoundaries = true,
		exclusionZones = {},
		mshFile = "examples/2D/damBreakKoshizuka/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "MinMax",
			outputFile = "tipPosition.txt",
			timeBetweenWriting = 0.01,
			minMax = "max",
			coordinate = 0 
		},
		{
			kind = "Mass",
			outputFile = "mass.txt",
			timeBetweenWriting = 0.01
		},
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.01,
			whatToWrite = {"p", "ke", "u", "v"},
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
		maxDT = 0.001,
		maxRemeshDT = -1,
		initialDT = 1e-8,
		securityCoeff = 0.1,
		
		MeshSmoother = {
			--[[{
				kind = "GETMe",
				epsTol = 0.01,
				maxIter = 20
			}]]--
		
		},
		
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
    local rhoStar = Problem.Material.rhoStar
    local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local g = -Problem.Solver.MomEq.bodyForce[2]
	local z0 = 2*0.146
	
	if(pos[2] <= z0 and pos[1] <= z0/2 + 1.1*Problem.Mesh.hchar) then
	    local rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(z0 - pos[2]) + 1)^(1/(K0p - 1))
		local p = K0/K0p*((rho/rhoStar)^K0p - 1)
		return {0, 0, p, rho, 0, 0}
	else 
		return {0, 0, 0, rhoStar, 0, 0}
	end
end

function Problem.Solver.MomEq.BC:BoundaryV(pos, t)
	return {0, 0}
end
