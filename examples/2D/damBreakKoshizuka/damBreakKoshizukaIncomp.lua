Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 4,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.0073,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.2,
		addOnFS = false,
		deleteFlyingNodes = true,
		laplacianSmoothingBoundaries = true, 
		boundingBox = {-0.01, -1, 0.594, 100},
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
		rho = 1000,
		gamma = 0
	},
	
	IC = {
		BoundaryFixed = true
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.001,
		maxRemeshDT = -1,
		initialDT = 0.001,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			gammaFS = 1,
			bodyForce = {0, -9.81},
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	local g = -Problem.Solver.MomContEq.bodyForce[2]
	local rho = Problem.Material.rho
	local p = 0
	if(pos[2] <= 2*0.146 and pos[1] <= 0.146 + 2*Problem.Mesh.hchar) then
		p = rho*g*(2*0.146 - pos[2])
	end
	return {0, 0, p}
end

function Problem.Solver.MomContEq.BC:BoundaryV(pos, t)
	return {0, 0}
end