Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 4,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.0048,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.5,
		addOnFS = false,
		deleteFlyingNodes = false,
		boundingBox = {-0.01, -0.01, 0.628, 100},
		exclusionZones = {{0.292, 0, 0.316, 0.048}},
		laplacianSmoothingBoundaries = false,
		mshFile = "examples/2D/damBreakWithObstacle/geometry.msh"
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
			timeBetweenWriting = 0.01,
		},
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
		maxDT = 0.0025,
		initialDT = 0.0025,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, -9.81},
			gammaFS = 0.5,
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	local rho = Problem.Material.rho
	local g = -Problem.Solver.MomContEq.bodyForce[2]
	local p
	if(pos[2] <= 2*0.146 and pos[1] <= 0.146 + 1.1*Problem.Mesh.hchar) then
		p = rho*g*(2*0.146 - pos[2])
	else 
		p = 0
	end
	return {0, 0, p}	
end

function Problem.Solver.MomContEq.BC:BoundaryV(pos, t)
	return {0, 0}
end