.PHONY: all clean phoniest

phoniest:
	poetry run seq2ids \
		--secrets_file local/secrets.yaml