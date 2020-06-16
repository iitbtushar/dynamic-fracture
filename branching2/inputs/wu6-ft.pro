log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};

control = 
{
  pause = 0.;
  runWhile = "i < 0";
};
  
input =
{
  file = "mesh.data";
};

// define the tolerance and max number of iteration
// for the AM solver (implemented in the ForkModule)
// AM_max_iteration = 1 => just one iteration as in Fenics code of Jian

amExp = 
{
  AM_tolerance     = 1e-4;
  AM_max_iteration = 1;
  convergence      = "Damage";
};


chain1 = 
{
  // note that: for this problem, this is just an empty file
  // as we do not need Dirichlet BCs for the phase field
  // you cannot remove the input so we need just an empty file
  input =
  {
    file = "phi.data";
  };

  model = "MP" // for parallel calculations
  {
    model = "Matrix"
    {
      matrix.type = "FEM";
      matrix.symmetric = true;
    
      model =  "VariationalDamPhase"
      {
          elements = "bulk";
          elemType = "Lagrange";
        
          shape =
          {
            type      = "Triangle3";
            intScheme = "Gauss1";
          };
    
          geom  = "Wu";

          material  = ["mat"];

          mat = 
          {
            elements = "bulk";
            degradation=  "WuDegradation"
            {
               softening="linear";
               young = 72e3;
               ft    = 28.65;
               Gc    = 3.8e-3;
               l0  = 0.25;

            };
          };
          
          writeCracks=
          {
            file     = "$(CASE_NAME)";
            interval = 1;
            x0       = 50.;
            y0       = 20.;
            phiMax   = 0.85;
          };
      };
    };
  };
  
  solver = "Nonlin"
  {
      precision  = 1.0e-5;
      //lineSearch = true;
      bounds     = ["b1"];

      solver = 
      {
          precision = 1.0e-7;
          type="CG";
          precon.type="ILUd";
          useThreads=true;
      };
      
      b1 = 
      {
         dofType    = "pf";
         lowerBound = 0.;
         upperBound = 1.;
      };
   };
  
  extraModules =
  {
    //modules =["energies", "view"];
    // use the following if you do not want to see the real time visualiation
    modules =["energies"];
    
    energies = "Sample"
    {
      file       = "$(CASE_NAME)-energies.dat";
      sampleWhen = "i % 2"; // every 4 steps
      header     = "Time | surface energy | strain energy | kinetic energy";
      dataSets   = [ "t","surfaceEnergy","strainEnergy", "kineticEnergy"];
    };
    

    view = "FemView"
    {
      window =
      {
        height = 300;
        width  = 900;
      };
    
      snapFile = "$(CASE_NAME)%2i.png";
      configFile  = "$(CASE_NAME).view";
    
      dataSets = [ "disp" ];
    
      disp =
      {
        type   = "Vector";
        vector = "state";
      };
    
      mesh =
      {
        //showFaces = false;
    
        plugins = "colors";
        colors =
        {
          type = "MeshColorView";
          data = "disp";
          palette   = "custom";
          autoScale = true;
        };
      };
    };
  };
};

chain2 = 
{
  input =
  {
    file = "u.data";
  };

  model = "MP"
  {
    model = "Matrix"
    {
      matrix.type = "FEM";
      matrix.symmetric = true;
    
      model       = "Multi"
      {
        models = [ "bulk", "top", "bot", "lodi" ];

        lodi =  "Lodi"
        {
          group = "top";
        };
    
        bulk = "VariationalDamMech"
        {
          elements = "bulk";
          elemType = "Lagrange";
    
          shape =
          {
            type      = "Triangle3";
            intScheme = "Gauss1";
          };
        
          density   = 2440e-12;
          E = 72e3; //to compute CFL
	      damThreshold = 1e-5; //for ref
    
          material = ["mat1"]; 
          mat1 = 
          {
            elements = "bulk";
            type   = "IsoPhaseField";
            dim    = 2;
            state  = "PLANE_STRESS";
    
            E     = 72e3;
            nu    = 0.22;
	        thickness = 1.0;
            
            degradation=  "WuDegradation"
            {
               softening="linear";
               young = 72e3;
               ft    = 28.65;
               Gc    = 3.8e-3;
               l0  = 0.25;
            };
            
            drivingForce = 
            {
              type = "Rankine";
            };
          };
        };
        
        top = 
        {
          type = "LoadScale";
    
          model =
          {
            type     = "LineLoad";
            elements = "notch-top";
            load     = -6.0;
            shape.type = "BLine2";
          };
        };

        bot = 
        {
          type = "LoadScale";
    
          model =
          {
            type     = "LineLoad";
            elements = "notch-bottom";
            load     = -6.0;
            shape.type = "BLine2";
          };
        };
      };
    };
  }; // end of MP model
  
  // Explicit dynamics
  
  solver=  
  {
    deltaTime = 3e-9; 
  };

  extraModules = 
  {
    modules = ["vtk", "sample"];

    vtk = "VTKOutput"
    {
        prefix   = "vtu/$(CASE_NAME)";
        format   = "binary"; // binary and appended 
        elements = "bulk";
        interval = 1000;
        data     = "stress|damage";
    };

    sample = "Sample"
    {
      file       = "$(CASE_NAME).dat";
      header     = "   t |  uy[mm] | fy[N] ";
      dataSets   = ["t", " chain2.model.model.model.lodi.disp[1]"," chain2.model.model.model.lodi.load[1]"]; 
    };

  };
};



