# Project-2-G3
ORF finder bacteria
while True:
        user_input = input("Enter your DNA sequence (spaces allowed): ")
    
    # Run your function
        result = validateSeq(user_input)

    if result is False:
        print("Invalid input! Only A, T, G and C are allowed.")
        print("-" * 60)
        continue

    final_seq = result
    break


print(f"\nYour validated sequence is: {final_seq}")

# for later