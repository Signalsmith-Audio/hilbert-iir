main: out/main
	@cd out && ./main
	
out/main: *.cpp ../*.h
	mkdir -p out
	g++ -std=c++11 -Ofast \
		-Wall -Wextra -Wfatal-errors -Wpedantic -pedantic-errors \
		main.cpp -o out/main

clean:
	rm -rf out
