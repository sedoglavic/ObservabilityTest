MAPLE=maple

release/maple.lib: src/*
	@echo "Cleaning the release folder"
	- @rm release/*
	@echo "Creating a new maple archive"
	@(echo "march('create',cat(currentdir(),\"/release/maple.lib\"),20);" | ${MAPLE} -q) 
	$(MAPLE) -q -i src/BaurStrassen.mpl
	$(MAPLE) -q -i src/dagnormal.mpl
	$(MAPLE) -q -i src/ObservabilityTest.mpl
	$(MAPLE) -q -i src/ObservabilitySymmetries.mpl

