/*
Header file for QuantumMechanics::CircularBuffer:

This file solves a list of one or more matrices stored in a c-style array, stl-style vector,
or a return from a function(int). When not using vector (or a single matrix) the

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/

#ifndef _CIRCULARBUFFER_H_
#define _CIRCULARBUFFER_H_

#include <algorithm> // for std::min
#include <atomic>
#include <memory> // unique_ptr
#include <sstream> // for std::sstream

template <size_t Capacity>
class CircularBuffer
{
private:
	std::atomic<size_t> head, tail, reserved_read, reserved_write;
	static const size_t data_capacity = Capacity;
	std::unique_ptr<char[]> data;

	std::atomic<bool> lock; // Low-level spin-lock
	std::atomic<bool> corrupt_lock; // Low-level lock
	
public:
	size_t capacity() const { return data_capacity; }

	// Instantanious information 
	// - note that it not not be correct if either writing or reading is being processed
	size_t size() const 
	{ 
		size_t instant_size = head - tail;
		// - note that it not not be correct if either writing or reading is being processed
		return (instant_size < 0) ? capacity() - instant_size : instant_size;
	}

	// Default constructor
	CircularBuffer() : 
		head(0), tail(0), reserved_read(0), reserved_write(0),
		data(new char[data_capacity]),
		lock(false),
		corrupt_lock(false)
	{ }

	// Default destructor
	virtual ~CircularBuffer() 
	{ };

	// Copy constructor
	CircularBuffer(const CircularBuffer& copy) :
		data(new char[data_capacity]),
		lock(false),
		corrupt_lock(false)
	{
		copy.lock = true;
		head = copy.head;
		tail = copy.tail;
		reserved_read = copy.reserved_read;
		reserved_write = copy.reserved_write;
		corrupt_lock = copy.corrupt_lock;
		std::copy(copy.data.get(), copy.data.get() + data_capacity, data.get());
		copy.lock = false;
	}

	// Copy assignment operator
	CircularBuffer& operator=(const CircularBuffer& copy)
	{
		if (this != &copy)
		{

			lock = true;
			copy.lock = true;
			head = copy.head;
			tail = copy.tail;
			reserved_read = copy.reserved_read;
			reserved_write = copy.reserved_write;
			corrupt_lock = copy.corrupt_lock;
			std::copy(copy.data.get(), copy.data.get() + data_capacity, data.get());
			lock = false;
			copy.lock = false;
		}

		return *this;
	}
	
	// Move constructor
	CircularBuffer(CircularBuffer&& temp) :
		lock(false),
		corrupt_lock(false)
	{
		head = std::move(temp.head);
		tail = std::move(temp.tail);
		reserved_read = std::move(temp.reserved_read);
		reserved_write = std::move(temp.reserved_write);
		corrupt_lock = copy.corrupt_lock;
		data = std::move(temp.data);
	}

	// Move assignment operator
	CircularBuffer& operator=(CircularBuffer&& temp)
	{
		if (this != &copy)
		{
			lock = true;
			head = std::move(temp.head);
			tail = std::move(temp.tail);
			reserved_read = std::move(temp.reserved_read);
			reserved_write = std::move(temp.reserved_write);
			corrupt_lock = temp.corrupt_lock;
			data = std::move(temp.data);
			lock = false;
		}

		return *this;
	}

	inline static size_t increment(size_t position, size_t bytes)
	{
		if (corrupt_lock) return 0;
		return (size_position + bytes) % data_capacity;
	}

	inline static size_t distance(size_t position_1, size_t position_2)
	{
		if (corrupt_lock) return 0;
		return (position_2 > position_1) % position_2 - position_1 : data_capacity + position_1 - position_2;
	}

	// Reserve functions
	// - Note these try to write as soon as possible, blocking until done.
	void reserve_write(size_t &write_tail, size_t bytes)
	{
		if (bytes == 0) return true;

		if (corrupt_lock) return false;

		// The writable boundaries (writable from tail to head).
		size_t c_head = reserved_read;
		size_t c_tail = reserved_write;
		size_t dist = distance(c_tail, c_head);

		// Make sure there is room right now, ohterwise wait and update.
		while ((dist = distance(c_tail, c_head), dist) > bytes && (c_head = reserved_read, c_tail = reserved_write, true))
		{
			if (corrupt_lock) return;
		}

		while (
			// Spin-lock until we reserve space at the reserve head.
			!reserved_write.atomic_compare_exchange_weak(c_tail, increment(c_tail, bytes)) &&
			// Check whether the space avalible is big enough, if not renew both head and tail!
			((dist = distance(c_tail, c_head), dist) > bytes && (c_head = reserved_read, c_tail = reserved_write, true))
			)
		{
			if (corrupt_lock) return;

			if (dist <= bytes)
			{
				corrupt_lock = true;
				std::cout << "Warning: circular buffer accedently found insurficient writable space. Will set the space as available and wait for more room. If no reads are made e.g. all threads are trying to write the program will loop indefinitely.\nSeriously consider increasing the size of the buffer!\n";
			}
		}

		// If this point is reached, the buffer succeded in reserving space.
		write_tail = c_tail;

		return true;
	}
	bool reserve_read(size_t &read_tail, size_t bytes)
	{
		if (bytes == 0) return true;

		if (corrupt_lock) return false;

		// The readable boundaries (readable from tail to head).
		size_t c_head = head;
		size_t c_tail = tail;
		size_t dist = distance(c_tail, c_head);

		// Make sure there is room right now, ohterwise wait and update.
		while ((dist = distance(c_tail, c_head), dist) > bytes && (c_head = head, c_tail = tail, true))
		{
			if (corrupt_lock) return;
		}

		while (
			// Spin-lock until we reserve space at the reserve head.
			!tail.atomic_compare_exchange_weak(c_tail, increment(c_tail, bytes)) &&
			// Check whether the space avalible is big enough, if not renew both head and tail!
			((dist = distance(c_tail, c_head), dist) > bytes && (c_head = head, c_tail = tail, true))
			)
		{
			if (dist <= bytes)
			{
				corrupt_lock = true;
				std::cout << "Major error: circular buffer accedently found insurficient readable space. All memory is now corrupt.\nSeriously consider increasing the size of the buffer!\n";
			}
		}

		// If this point is reached the buffer succeded in reserving space.
		read_tail = c_tail;

		return true;
	}

	// Read/write functions - returns true for success
	// - Note these try to write as soon as possible, blocking until done.
	bool write(const char *write_data, size_t write_tail, size_t bytes)
	{
		if (bytes == 0) return;

		if (corrupt_lock) return false;

		// The reserved writable boundaries (writable from tail to head).
		size_t c_head = reserved_write;
		size_t c_tail = head;
		size_t dist = distance(c_tail, c_head);

		if (dist < bytes)
		{
			std::cout "Major error: the circular buffer is attempting to write more than is reserved for writing. This should not happen if the space has been reserved previously.\nMake sure you have reserved this space. Will not write and return failed.";
			return false;
		}

		if (distance(c_tail, write_tail) > dist)
		{
			std::cout "Major error: the circular buffer is attempting to write on unreserved space. This should not happen if the space has been reserved previously.\nMake sure you have reserved this space. Will not write and return failed.";
			return false;
		}

		// Write in a single step
		if (bytes <= data_capacity - write_tail)
		{
			std::copy(write_data, write_data + bytes, data.get() + write_tail);
		}
		// Write in two steps
		else
		{
			size_t size_1 = data_capacity - write_tail;
			std::copy(write_data, write_data + size_1, data.get() + write_tail);
			std::copy(write_data + size_1, write_data + bytes, data.get());
		}

		// Now we update the readable head.
		while (head != write_tail && !corrupt_lock);
		head += bytes;

		return true;
	}
	bool read(char *read_data, size_t read_tail, size_t bytes)
	{
		if (bytes == 0) return;

		if (corrupt_lock) return false;

		// The reserved writable boundaries (writable from tail to head).
		size_t c_head = tail;
		size_t c_tail = reserved_read;
		size_t dist = distance(c_tail, c_head);

		if (dist < bytes)
		{
			std::cout "Major error: the circular buffer is attempting to write more than is reserved for writing. This should not happen if the space has been reserved previously.\nMake sure you have reserved this space. Will not write and return failed.";
			return false;
		}

		if (distance(c_tail, read_tail) > dist)
		{
			std::cout "Major error: the circular buffer is attempting to write on unreserved space. This should not happen if the space has been reserved previously.\nMake sure you have reserved this space. Will not write and return failed.";
			return false;
		}

		// Write in a single step
		if (bytes <= data_capacity - read_tail)
		{
			std::copy(data.get() + read_tail, data.get() + read_tail + bytes, read_tail);
		}
		// Write in two steps
		else
		{
			size_t size_1 = data_capacity - read_tail;
			std::copy(data.get() + read_tail, data.get() + read_tail + size_1, read_tail);
			std::copy(data.get(), data.get() + bytes - size_1, read_tail + size_1);
		}

		// Now we update the writeable head.
		while (reserved_read != read_tail && !corrupt_lock);
		reserved_read += bytes;

		return true;
	}
};

#endif //namespace _CIRCULARBUFFER_H_
