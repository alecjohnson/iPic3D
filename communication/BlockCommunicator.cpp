#include "BlockCommunicator.h"
#include "Parameters.h"

BlockCommunicator::BlockCommunicator(Connection connection_)
{
  commState=NONE;
  init(connection_);
}

void BlockCommunicator::init(Connection connection_)
{
  assert(commState==NONE);
  connection = connection_;
  commState = INITIAL;
  blocksize = Parameters::get_blockSize();
  numblocks = Parameters::get_numBlocks();

  assert(numblocks>0);
  for(nextid=0;nextid<numblocks;nextid++)
  {
    Block<type>* newBlock = new Block<type>(blocksize, nextid);
    blockList.push_back(newBlock);
  }
  curr_block = blockList.begin();
}

BlockCommunicator::~BlockCommunicator()
{
  std::list<void*>::iterator i;
  for(i=blockList.begin(); i != blockList.end(); ++i)
    delete &fetch_block(i);
}
