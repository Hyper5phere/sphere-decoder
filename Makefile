PROGRAM=sphdec

subsystem:
	$(MAKE) -C src

.PHONY: clean
clean:
	cd src && $(MAKE) clean

run:
	./$(PROGRAM)

debug:
	gdb $(PROGRAM)
