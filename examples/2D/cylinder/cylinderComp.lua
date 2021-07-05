Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 15,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.01,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.55,
		addOnFS = true,
		deleteFlyingNodes = false,
		boundingBox = {-0.1, -0.52, 1.05, 0.52},
		exclusionZones = {},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/2D/cylinder/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.5,
			whatToWrite = {"u", "v", "p", "ke", "velocity"},
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
		gamma = 0,
		K0 = 22000,
		K0p = 7.6,
		rhoStar = 1000
	},
	
	IC = {
		WallsFixed = true,
		CylinderWallFixed = true,
		FluidInputFixed = true
	},
	
	Solver = {
	    id = "CDS_dpdt",
		adaptDT = true,
		maxDT = 0.05,
		maxRemeshDT = -1,
		initialDT = 1e-8,
		securityCoeff = 0.1,
		
		MeshSmoother = {
			--[[{
				kind = "PSmoother",
				a = 0.001,
				epsADRTol = 0.001,
				betaInit = 10
			},
			{
				kind = "GETMe",
				epsTol = 0.01,
				maxIter = 20
			}--]]
		},
		
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
	return {0, 0, 0, rhoStar, 0, 0}
end

function Problem.IC:initFluidInputStates(pos)
	local rhoStar = Problem.Material.rhoStar
	return {0, 0, 0, rhoStar, 0, 0}
end

function Problem.IC:initCylinderWallStates(pos)
	local rhoStar = Problem.Material.rhoStar
	return {0, 0, 0, rhoStar, 0, 0}
end

function Problem.Solver.MomEq.BC:WallsV(pos, t)
	return {0, 0}
end

function Problem.Solver.MomEq.BC:CylinderWallV(pos, t)
	return {0, 0}
end

function Problem.Solver.MomEq.BC:FluidInputV(pos, t)
	local DT = 2
	if(t > DT) then
		return {0, 0}
	else
		local v0 = 0.01*4
		local acc = -6*v0*t^2/DT^3 + 6*v0*t/DT^2
		return {acc*(0.375^4 - pos[2]^4)/0.375^4, 0}
	end
end