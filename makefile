SRCFILES=$(shell find -maxdepth 1 -name '*.c')
TARGET= $(basename $(SRCFILES))

all : $(TARGET)

define outputrule
$(1): $(1:%=%.c)
	gcc -o $$@ -Wall -O2 -pipe -march=native $$< -lm -lpthread
endef
$(foreach src,$(TARGET),$(eval $(call outputrule,$(src))))
