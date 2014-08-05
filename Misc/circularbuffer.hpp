/*
Header file for QuantumMechanics::CircularBuffer:

This file declares and defines the circular buffer (ring buffer).
The particular implementation is mostly! lock-free multi-producer-multi-consumer atomic buffer.

By mostly, we mean that reads (writes) try to lock in a reserved space, which is not nessecarily immediate.

---
Copyright (C) 2014, Søren Schou Gregersen <sorge@nanotech.dtu.dk>
*/

#ifndef _CIRCULARBUFFER_H_
#define _CIRCULARBUFFER_H_

#include <algorithm> // for std::min
#include <atomic> // for std::atomic
#include <memory> // for std::unique_ptr
#include <iostream> // for std::cout

template <size_t Capacity>
class CircularBuffer
{
private:
	// enum to declare space available for read/write or either just written or just read.
	// - This is needed to update newly written/read sectors that can't immidiately be declared available.
	enum class SpaceState {
		Neutral,
		RecentlyWritten,
		RecentlyRead
	};

	std::atomic<size_t> head, tail, reserved_read, reserved_write; // position identifiers
	static const size_t data_capacity = Capacity; // Constant size
	std::unique_ptr<unsigned char[]> data; // The storage space
	std::unique_ptr<SpaceState[]> spaces; // The declaration space for inbetween simultanous read/writes.
	
public:
	// Instantanious information 
	size_t capacity() const { return data_capacity; }

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
		spaces(new SpaceState[data_capacity])
	{
		for (auto& space : (SpaceState(*)[data_capacity]) spaces)
			space = SpaceState::Neutral;
	}

	// Default destructor
	virtual ~CircularBuffer() 
	{ };

	// Copy constructor
	// - note that lock-free means no protection against simultanious reading or writing.
	CircularBuffer(const CircularBuffer& copy) :
		data(new char[data_capacity]),
		spaces(new SpaceState[data_capacity])
	{
		head = copy.head;
		tail = copy.tail;
		reserved_read = copy.reserved_read;
		reserved_write = copy.reserved_write;
		std::copy(copy.data.get(), copy.data.get() + data_capacity, data.get());
		std::copy(copy.spaces.get(), copy.spaces.get() + data_capacity, spaces.get());
	}

	// Copy assignment operator
	// - note that lock-free means no protection against simultanious reading or writing.
	CircularBuffer& operator=(const CircularBuffer& copy)
	{
		if (this != &copy)
		{
			head = copy.head;
			tail = copy.tail;
			reserved_read = copy.reserved_read;
			reserved_write = copy.reserved_write;
			std::copy(copy.data.get(), copy.data.get() + data_capacity, data.get());
			std::copy(copy.spaces.get(), copy.spaces.get() + data_capacity, spaces.get());
		}

		return *this;
	}
	
	// Move constructor
	CircularBuffer(CircularBuffer&& temp)
	{
		head = std::move(temp.head);
		tail = std::move(temp.tail);
		reserved_read = std::move(temp.reserved_read);
		reserved_write = std::move(temp.reserved_write);
		data = std::move(temp.data);
		spaces = std::move(temp.spaces);
	}

	// Move assignment operator
	CircularBuffer& operator=(CircularBuffer&& temp)
	{
		if (this != &copy)
		{
			head = std::move(temp.head);
			tail = std::move(temp.tail);
			reserved_read = std::move(temp.reserved_read);
			reserved_write = std::move(temp.reserved_write);
			data = std::move(temp.data);
			spaces = std::move(temp.spaces);
		}

		return *this;
	}

protected:
	// Helper function
	// - find the position accounting for ring buffer wrap.
	inline static size_t increment(size_t position, size_t bytes)
	{
		return (size_position + bytes) % data_capacity;
	}

	// Helper function
	// - find the distance accounting for ring buffer wrap (along positive direction).
	inline static size_t distance(size_t position_1, size_t position_2)
	{
		return (position_2 > position_1) % position_2 - position_1 : data_capacity + position_1 - position_2;
	}

public:
	// Read/write functions
	// - Note these try to reserve space, atomically determining if there is space and reserving.
	// - If the space becomes reserved the read/write happens afterwards. Unfortunately, the function blocks until it can unreserve the space again. Returns true on success.
	// - If there is no space found, it return false.
	bool write(const unsigned char *write_data, size_t bytes)
	{
		if (bytes == 0) return true; // No need.

		// The writable boundaries (writable from tail to head).
		size_t c_head = reserved_read;
		size_t c_tail = reserved_write;

		// Make sure there is room right now, ohterwise return failed.
		do
		{
			c_head = reserved_read; // update the head (in case there is more room).
			size_t dist = distance(c_tail, c_head);
			// note the head might have changed already, leaving more room. We disregard this situation.

			if (dist < bytes) // No room for the data.
			{
				std::cout << "A circular buffer has become full.\n As this slows everything, seriously consider making the capacity larger!" << std::endl;
				return false;
			}

		} while (!reserved_write.atomic_compare_exchange_weak(c_tail, increment(c_tail, bytes)));
		// if the tail can be locked, it is incremented. Note the distance is already tested and approved.
		// This is the only spin-lock used, and it should not last that long.

		// If this point is reached, the buffer succeded in reserving space.

		// Write in a single step
		if (bytes <= data_capacity - c_tail)
		{
			std::copy(write_data, write_data + bytes, data.get() + c_tail);

			for (size_t i = c_tail; i < c_tail + bytes; i++)
				spaces[i] = SpaceState::RecentlyWritten;
		}
		// Write in two steps
		else
		{
			size_t size_1 = data_capacity - c_tail;
			std::copy(write_data, write_data + size_1, data.get() + c_tail);
			for (size_t i = c_tail; i < c_tail + size_1; i++)
				spaces[i] = SpaceState::RecentlyWritten;
			std::copy(write_data + size_1, write_data + bytes, data.get());
			for (size_t i = 0; i < bytes - size_1; i++)
				spaces[i] = SpaceState::RecentlyWritten;
		}

		// Now, if this write is not the first write after head, we do nothing.
		if (head != c_tail) return true;

		// Otherwise, we would like to declare this space ready for read.
		// - Note that all other treads exit as they are not at the head, so it will remain c_tail until this thread changes it.

		size_t new_head = head;

		while (spaces[new_head] == SpaceState::RecentlyWritten) // For each written space (also from other writes)...
		{
			spaces[new_head] = SpaceState::Neutral; // Declare it neutral (not recently written or read).
			new_head = increment(new_head, 1); // Increment the new head.
		}

		// Now we update the readable head.
		head = new_head;

		return true;
	}
	bool read(unsigned char *read_data, size_t bytes)
	{
		if (bytes == 0) return true; // No need.

		// The readable boundaries (readable from tail to head).
		size_t c_head = head;
		size_t c_tail = tail;

		// Make sure there is data right now, ohterwise return failed.
		do
		{
			c_head = head; // update the head (in case there is more data now).
			size_t dist = distance(c_tail, c_head);
			// note the head might have changed already, leaving more data. We disregard this situation.

			if (dist < bytes) // Not enough data.
			{
				std::cout << "A circular buffer has been asked to supply more data than it can at the moment.\n This may be on purpose, however this is usually not the case." << std::endl;
				return false;
			}

		} while (!tail.atomic_compare_exchange_weak(c_tail, increment(c_tail, bytes))); 
		// if the tail can be locked, it is incremented. Note the distance is already tested and approved.
		// This is the only spin-lock used, and it should not last that long.
			
		// If this point is reached the buffer succeded in reserving space.

		// Write in a single step
		if (bytes <= data_capacity - c_tail)
		{
			std::copy(data.get() + c_tail, data.get() + c_tail + bytes, read_data);
			for (size_t i = c_tail; i < c_tail + bytes; i++)
				spaces[i] = SpaceState::RecentlyRead;
		}
		// Write in two steps
		else
		{
			size_t size_1 = data_capacity - c_tail;
			std::copy(data.get() + c_tail, data.get() + c_tail + size_1, read_data);
			for (size_t i = c_tail; i < c_tail + size_1; i++)
				spaces[i] = SpaceState::RecentlyRead;
			std::copy(data.get(), data.get() + bytes - size_1, read_data + size_1);
			for (size_t i = 0; i < bytes - size_1; i++)
				spaces[i] = SpaceState::RecentlyRead;
		}

		// Now, if this read is not the first read after free space, we do nothing.
		if (reserved_read != c_tail) return true;

		// Otherwise, we would like to declare this space ready for write.
		// - Note that all other treads exit as they are not at the head, so it will remain c_tail until this thread changes it.

		size_t new_head = reserved_read;

		while (spaces[new_head] == SpaceState::RecentlyRead) // For each read space (also from other reads)...
		{
			spaces[new_head] = SpaceState::Neutral; // Declare it neutral (not recently written or read).
			new_head = increment(new_head, 1); // Increment the new head.
		}

		// Now we update the writeable head.
		reserved_read = new_head;

		return true;
	}
};

#endif //namespace _CIRCULARBUFFER_H_
